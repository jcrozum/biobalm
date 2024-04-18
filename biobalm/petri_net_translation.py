"""
Utilities for translating Boolean networks into Petri nets.

Implements the translation from a `BooleanNetwork` object into a Petri net that
can be processed by Trappist (for finding trap spaces). The Petri net is
represented as a `networkx.DiGraph`, with nodes having either a `kind=place` or
a `kind=transition` attribute.

The variable names in the network have to be "sanitized" before translation. In
particular, this means they can't use any special characters beyond "_". In the
resulting Petri net, we then use a "b0_*" and "b1_*" prefix to distinguish the
"zero" and "one" places created for each variable. This is also important for
`clingo`, as we have to guarantee that in our logic program, symbols start with
a lowercase letter.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Generator

    from biodivine_aeon import Bdd, BooleanNetwork

    from biobalm.types import BooleanSpace

import copy
import re
from typing import cast

from biodivine_aeon import (
    BddPartialValuation,
    BddVariableSet,
    BoolType,
    SymbolicContext,
)
from networkx import DiGraph  # type: ignore

# Enables statistics logging.
DEBUG = False
"""Enables debug logging to stdout."""


def sanitize_network_names(network: BooleanNetwork, check_only: bool = False):
    """
    Rename variables in a network so that they can be safely used in Trappist.

    Verifies that all network variable names contain only alphanumeric
    characters and underscores. If this is not the case, attempts to rename the
    variables to make them compliant. Returns a *copy* of the original network
    that only uses sanitized names.

    Note that AEON should already prune away most of the special characters when
    parsing models, but it still allows names like `Gene_{Subscript}` (i.e.
    curly brackets) which we have to sanitize here.

    If `check_only=True` is specified, no renaming takes place and the function
    fails with a `RuntimeError` instead.

    Parameters
    ----------
    network : BooleanNetwork
        The network to sanitize.
    check_only : bool
        If `True`, no renaming takes place and the function fails with a
        `RuntimeError` if the network contains invalid variable names.

    Returns
    -------
    BooleanNetwork
        A copy of the original network with sanitized variable names.
    """
    network = copy.copy(network)
    for var in network.variables():
        name = network.get_variable_name(var)
        if not re.match("^[a-zA-Z0-9_]+$", name):
            if check_only:
                raise RuntimeError(f"Found unsanitized variable: `{name}`.")
            # Replace all invalid characters with an underscore
            new_name = re.sub("[^a-zA-Z0-9_]", "_", name)
            while True:
                try:
                    network.set_variable_name(var, new_name)
                    break
                except Exception:
                    # Name clash happened. Try to resolve this by adding an extra underscore to
                    # the variable name. In theory, this can repeat until a unique name is found.
                    new_name = "_" + new_name

    return network


def variable_to_place(variable: str, positive: bool) -> str:
    """
    Generate a Petri net place name from a network variable name.

    Parameters
    ----------
    variable : str
        The name of the network variable.
    positive : bool
        `True` if the place corresponding to the variable should be positive,
        `False` if it shoudl be negative.

    Returns
    -------
    str
        The name of the corresponding Petri net place.
    """
    if positive:
        return f"b1_{variable}"
    else:
        return f"b0_{variable}"


def place_to_variable(place: str) -> tuple[str, bool]:
    """
    Extract the variable name and state from a Petri net place name.

    Parameters
    ----------
    place : str
        The name of the Petri net place.

    Returns
    -------
    tuple[str, bool]
        The name of the variable, and whether the variable is positive or negative.
    """
    if place.startswith("b1_"):
        return (place[3:], True)
    elif place.startswith("b0_"):
        return (place[3:], False)
    else:
        raise Exception(f"Invalid place name: `{place}`.")


def extract_variable_names(encoded_network: DiGraph) -> list[str]:
    """
    Extract the variable names from a Petri net encoded Boolean network.

    The variables are  sorted lexicographically, since the original BN ordering
    is not preserved by the Petri net. However, BNs order variables
    lexicographically by default, so unless the Petri net was created from a
    custom BN (i.e. not from a normal model file), the ordering should be the
    same.

    Parameters
    ----------
    encoded_network : DiGraph
        The Petri net encoded Boolean network.

    Returns
    -------
    list[str]
        The list of variable names.
    """
    variables: list[str] = []
    for node in encoded_network.nodes():  # type: ignore
        node = str(node)  # type: ignore
        if node.startswith("b0_"):
            variables.append(place_to_variable(node)[0])

    return sorted(variables)


def extract_source_variables(encoded_network: DiGraph) -> list[str]:
    """
    List variable names that represent source nodes of the encoded network.

    Source nodes are those nodes with an identity update function.

    Parameters
    ----------
    encoded_network : DiGraph
        The Petri net encoded Boolean network.

    Returns
    -------
    list[str]
        The list of source variable names.
    """
    variables = extract_variable_names(encoded_network)
    source_set = set(variables)
    for _, change_var in encoded_network.nodes(data="change"):  # type: ignore
        if change_var in source_set:
            source_set.remove(change_var)  # type: ignore[reportUnknownArgumentType] # noqa
    source_nodes: list[str] = sorted(source_set)
    return source_nodes


def restrict_petrinet_to_subspace(
    petri_net: DiGraph,
    sub_space: BooleanSpace,
) -> DiGraph:
    """
    Create a copy of a Petri net restricted to a sub-space.

    Creates a copy of the given Petri net, but with the variables given in
    `sub_space` fixed to their respective values.

    Note that this completely eliminates the constant variables from the Petri
    net, but it does not perform any further constant propagation or
    percolation. Variables that are fixed in the `sub_space` but do not exist in
    the Petri net are ignored.

    The `sub_space` can contain variables that do not appear
    in the `petri_net`. Such variables are simply ignored.

    Parameters
    ----------
    petri_net : DiGraph
        The Petri net to restrict.
    sub_space : BooleanSpace
        The sub-space to restrict the Petri net to.

    Returns
    -------
    DiGraph
        The restricted Petri net.
    """
    result = copy.deepcopy(petri_net)
    for var, value in sub_space.items():
        # The idea is that we want to *remove* the place that corresponds to the fixed
        # value (it's effect on transitions is assumed to be fulfilled). Then, we remove
        # all transitions that depend on the inverse of the fixed value, as these can
        # never be satisfied. Then, we also remove inverse place.
        fixed_place = variable_to_place(var, bool(value))
        inverse_place = variable_to_place(var, not bool(value))

        if fixed_place not in result.nodes or inverse_place not in result.nodes:
            # These nodes were already removed.
            continue

        # First, remove all transitions that modify the fixed variable.
        # These are removed regardless of the actual value.

        to_delete: set[str] = set()

        # Transitions that put marker into one of the place values, but do not take it.
        for tr in result.predecessors(fixed_place):  # type: ignore
            if not result.has_edge(fixed_place, tr):  # type: ignore
                to_delete.add(cast(str, tr))
        for tr in result.predecessors(inverse_place):  # type: ignore
            if not result.has_edge(inverse_place, tr):  # type: ignore
                to_delete.add(cast(str, tr))

        # Transitions that depend on the value of the inverse place
        for tr in result.successors(inverse_place):  # type: ignore
            to_delete.add(cast(str, tr))

        for tr in to_delete:
            result.remove_node(tr)  # type: ignore

        result.remove_node(fixed_place)  # type: ignore
        result.remove_node(inverse_place)  # type: ignore
    return result


def network_to_petrinet(
    network: BooleanNetwork, symbolic_context: SymbolicContext | None = None
) -> DiGraph:
    """
    Convert a Boolean network to a Petri net.

    Converts a `biodivine_aeon.BooleanNetwork` to a `DiGraph` representing a Petri net encoding
    of the original network. For details about the encoding, see module
    description.

    Note that the given network needs to have "sanitized" names, otherwise the
    method will fail (see `sanitize_network_names` in this module).

    The operation uses translation through `biodivine_aeon.Bdd` to generate a
    disjunctive normal form of the network's update functions. This is
    facilitated by `biodivine_aeon.SymbolicContext`. If a suitable context already
    exists, it can be provided as the second argument. Otherwise it will be
    created.

    Parameters
    ----------
    network : BooleanNetwork
        The network to convert.
    symbolic_context : SymbolicContext | None
        The context used for the symbolic conversion, as an
        `biodivine_aeon.SymbolicContext` object. This is a mapping from the
        network nodes to BDD variables that preserves variable ordering in BDDs.
        If not given, a new one will be created from the network.

    Returns
    -------
    DiGraph
        The Petri net encoding of the given network.
    """
    # Assert that all network names are already sanitized.
    sanitize_network_names(network, check_only=True)

    assert (
        network.explicit_parameter_count() == 0
    ), f"Parametrized networks are not supported. Found parameters: {network.explicit_parameter_names()}."

    # Implicit parameters with no regulators are allowed, since they just reprtesent free inputs
    # and are explicitly handled by the succession diagram.
    non_input_implicit = [
        v for v in network.implicit_parameters() if len(network.predecessors(v)) > 0
    ]
    if len(non_input_implicit) > 0:
        names = [network.get_variable_name(x) for x in non_input_implicit]
        raise AssertionError(
            f"Parametrized networks are not supported. Found implicit parameters: {names}."
        )

    if symbolic_context is None:
        symbolic_context = SymbolicContext(network)

    pn = DiGraph()

    # Create a positive and a negative place for each variable.
    places: dict[str, tuple[str, str]] = {}
    for name in network.variable_names():
        p_name = variable_to_place(name, positive=True)
        n_name = variable_to_place(name, positive=False)
        pn.add_node(p_name, kind="place")  # type: ignore[reportUnknownMemberType]
        pn.add_node(n_name, kind="place")  # type: ignore[reportUnknownMemberType]
        places[name] = (n_name, p_name)

    # Create PN transitions for implicants of every BN transition.
    for var in network.variables():
        var_name = network.get_variable_name(var)
        update_function = network.get_update_function(var)
        if update_function is None:
            # This variable is a constant input with a free update function.
            assert len(network.predecessors(var)) == 0
            continue

        function_bdd = symbolic_context.mk_update_function(update_function)
        var_bdd = symbolic_context.mk_network_variable(var)

        p_bdd = function_bdd.l_and(var_bdd.l_not())
        n_bdd = function_bdd.l_not().l_and(var_bdd)

        if DEBUG:
            print(f"Start translation for `{var_name}`: {len(p_bdd)} | {len(n_bdd)}")

        # Add 0->1 edges.
        _create_transitions(
            pn, symbolic_context.bdd_variable_set(), places, var_name, p_bdd, go_up=True
        )
        # Add 1-> 0 edges.
        _create_transitions(
            pn,
            symbolic_context.bdd_variable_set(),
            places,
            var_name,
            n_bdd,
            go_up=False,
        )

    return pn


def optimized_recursive_dnf_generator(
    bdd: Bdd,
) -> Generator[BddPartialValuation, None, None]:
    """
    Yields a generator of `BddPartialValuation` objects, similar to
    `Bdd.clause_iterator` in AEON, but uses a recursive optimization strategy
    to return a smaller result than the default method `Bdd` clause sequence.
    Note that this is still not the "optimal" DNF, but is often close enough.

    This is technically slower for BDDs that already have a small clause count,
    but can be much better in the long-term when the clause count is high.
    """

    # Some performance notes:
    #  - We could cache the best restriced BDDs, but usually this does not help much.
    #  - We could try optimizing based on clause count instead of BDD size, but this will
    #    need a new method in `lib-bdd` to compute clause count efficiently. (Now
    #    we can only do it by walking the iterator.)
    #  - The initial BDD size/ordering still do matter: starting with a smaller but equivalent
    #    BDD can help achieve smaller DNF size.

    if bdd.is_false():
        return
    if bdd.is_true():
        yield BddPartialValuation(bdd.__ctx__(), cast(dict[str, BoolType], {}))
        return

    support = sorted(bdd.support_set())
    best_var = support[0]
    best_size = 10 * len(bdd)
    for var in support:
        t = bdd.r_restrict({var: True})
        f = bdd.r_restrict({var: False})
        t_size = len(t)
        f_size = len(f)
        if t_size + f_size < best_size:
            best_size = t_size + f_size
            best_var = var

    for t_val in optimized_recursive_dnf_generator(bdd.r_restrict({best_var: True})):
        t_val[best_var] = True
        yield t_val

    for f_val in optimized_recursive_dnf_generator(bdd.r_restrict({best_var: False})):
        f_val[best_var] = False
        yield f_val


def _create_transitions(
    pn: DiGraph,
    ctx: BddVariableSet,
    places: dict[str, tuple[str, str]],
    var_name: str,
    implicant_bdd: Bdd,
    go_up: bool,
):
    """
    Just a helper method that creates PN transitions from BDDs representing
    positive/negative BN transitions.
    """
    dir_str = "up" if go_up else "down"
    total = 0
    for t_id, implicant in enumerate(optimized_recursive_dnf_generator(implicant_bdd)):
        total += 1
        t_name = f"tr_{var_name}_{dir_str}_{t_id + 1}"
        pn.add_node(  # type: ignore[reportUnknownMemberType]
            t_name,
            kind="transition",  # type: ignore[reportUnknownMemberType]
            change=var_name,
            direction=dir_str,
        )
        # The transition moves a token either from "zero place" to the
        # "one place", or vice versa.
        if go_up:
            pn.add_edge(places[var_name][0], t_name)  # type: ignore[reportUnknownMemberType] # noqa
            pn.add_edge(t_name, places[var_name][1])  # type: ignore[reportUnknownMemberType] # noqa
        else:
            pn.add_edge(places[var_name][1], t_name)  # type: ignore[reportUnknownMemberType] # noqa
            pn.add_edge(t_name, places[var_name][0])  # type: ignore[reportUnknownMemberType] # noqa
        for variable, value in implicant.items():
            variable_str = ctx.get_variable_name(
                variable
            )  # Convert from BDD variable to name.
            if variable_str == var_name:
                continue
            # For the remaining variables, we simply check if the required
            # token is present in the corresponding place.
            pn.add_edge(places[variable_str][value], t_name)  # type: ignore[reportUnknownMemberType] # noqa
            pn.add_edge(t_name, places[variable_str][value])  # type: ignore[reportUnknownMemberType] # noqa
    if DEBUG:
        print(f"  >> Generated {total} total PN transitions.")
