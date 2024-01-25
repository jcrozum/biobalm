from __future__ import annotations

"""
    Implements the translation from a `BooleanNetwork` object into a Petri net that can be processed
    by Trappist. The Petri net is represented as a `DiGraph`, with nodes having either a `kind=place`
    or a `kind=transition` attribute.

    The variable names in the network have to be "sanitized" before translation. In particular, this
    means they can't use any special characters beyond "_". In the resulting Petri net, we then use
    a "b0_*" and "b1_*" prefix to distinguish the "zero" and "one" places created for each variable.
    This is also important for `clingo`, as we have to guarantee that in our logic program, symbols
    start with a lowercase letter.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from biodivine_aeon import BooleanNetwork, SymbolicContext, Bdd

import re
import copy

from networkx import DiGraph  # type: ignore


def sanitize_network_names(network: BooleanNetwork, check_only: bool = False):
    """
    Verifies that all network variable names contain only alphanumeric characters and "_".
    If this is not the case, attempts to rename the variables to make them compliant.
    Returns a *copy* of the original network that only uses sanitized names.

    Note that AEON should already prune away most of the special characters when parsing models,
    but it still allows names like `Gene_{Subscript}` (i.e. curly brackets) which we have
    to sanitize here.

    If `check_only=True` is specified, no renaming takes place and the function fails with
    a `RuntimeError` instead.
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
                except Exception:
                    # Name clash happened. Try to resolve this by adding an extra underscore to 
                    # the variable name. In theory, this can repeat until a unique name is found.
                    new_name = "_" + new_name                    

    return network


def variable_to_place(variable: str, positive: bool) -> str:
    """
    Convert the name of a network variable to the name of the corresponding positive
    or negative Petri net place.
    """
    if positive:
        return f"b1_{variable}"
    else:
        return f"b0_{variable}"


def place_to_variable(place: str) -> tuple[str, bool]:
    """
    Convert the name of a Petri net place to the name of the network variable, plus
    a Boolean indicating whether the original place was positive or negative.
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

    The variables are  sorted lexicographically, since the original BN ordering is not
    preserved by the Petri net. However, BNs order variables lexicographically by default,
    so unless the Petri net was created from a custom BN (i.e. not from a normal model file),
    the ordering should be the same.
    """
    variables: list[str] = []
    for node in encoded_network.nodes():  # type: ignore
        node = str(node)  # type: ignore
        if node.startswith("b0_"):
            variables.append(place_to_variable(node)[0])

    return sorted(variables)


def network_to_petrinet(network: BooleanNetwork, ctx: SymbolicContext | None = None) -> DiGraph:
    """
    Converts a `BooleanNetwork` to a `DiGraph` representing a Petri net encoding
    of the original network. For details about the encoding, see module description.

    Note that the given network needs to have "sanitized" names, otherwise the
    method will fail (see `sanitize_network_names` in this module).

    The operation uses translation through `biodivine_aeon.Bdd` to generate a disjunctive normal form 
    of the network's update functions. This is facilitated by `biodivine_aeon.SymbolicContext`. If such
    context already exists, it can be provided as the second argument. Otherwise it will be created.
    """
    # Assert that all network names are already sanitized.
    sanitize_network_names(network, check_only=True)

    assert (
        network.explicit_parameter_count() == 0
    ), "Unsupported: Network contains explicit parameters."
    assert (
        network.implicit_parameter_count() == 0
    ), "Unsupported: Network contains implicit parameters."

    if ctx is None:
        ctx = SymbolicContext(network)

    pn = DiGraph()

    # Create a positive and a negative place for each variable.
    places = {}
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
        # All function must be known if there are no parameters.
        assert update_function is not None 
        function_bdd = ctx.mk_update_function(update_function)
        var_bdd = ctx.mk_network_variable(var)

        positive_implicants = function_bdd.l_and(var_bdd.l_not())
        negative_implicants = function_bdd.l_not().l_and(var_bdd)

        # Add 0->1 edges.
        _create_transitions(pn, ctx, places, var_name, positive_implicants, go_up=True)
        # Add 1-> 0 edges.
        _create_transitions(pn, ctx, places, var_name, negative_implicants, go_up=False)

    return pn


def _create_transitions(
    pn: DiGraph,
    ctx: SymbolicContext,
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
    for t_id, implicant in enumerate(implicant_bdd.clause_iterator()):
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
            variable_str = ctx.get_network_variable_name(variable)  # Convert from BDD variable to name.
            if variable_str == var_name:
                continue
            # For the remaining variables, we simply check if the required
            # token is present in the corresponding place.
            pn.add_edge(places[variable_str][value], t_name)  # type: ignore[reportUnknownMemberType] # noqa
            pn.add_edge(t_name, places[variable_str][value])  # type: ignore[reportUnknownMemberType] # noqa
