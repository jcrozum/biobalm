from __future__ import annotations

"""
    Here, we implement the Trappist method for computing fixed-points, minimum trap spaces
    and maximum trap spaces, including time-reversed networks.
"""
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Callable
    from clingo import Model

from biodivine_aeon import BooleanNetwork
from clingo import Control, SolveHandle
from networkx import DiGraph  # type: ignore

from nfvsmotifs.petri_net_translation import (
    extract_variable_names,
    network_to_petrinet,
    place_to_variable,
    variable_to_place,
)


def trappist_async(
    network: BooleanNetwork | DiGraph,
    on_solution: Callable[[dict[str, int]], bool],
    problem: str = "min",
    reverse_time: bool = False,
    ensure_subspace: dict[str, int] | None = None,
    avoid_subspaces: list[dict[str, int]] | None = None,
):
    """
    The same as the `trappist` method, but instead of returning a list of spaces as a result, the
    spaces are returned to the supplied `on_solution` callback. You can stop the enumeration by
    returning `False` from this callback.
    """
    if ensure_subspace is None:
        ensure_subspace = {}
    if avoid_subspaces is None:
        avoid_subspaces = []

    if isinstance(network, BooleanNetwork):
        bn = network
        petri_net = network_to_petrinet(network)
    else:
        bn = None
        petri_net = network

    assert isinstance(petri_net, DiGraph)

    if bn is None:
        variables = extract_variable_names(petri_net)
    else:
        variables = [bn.get_variable_name(v) for v in bn.variables()]

    # Source node is a node that has no transitions in the PN encoding
    # (i.e. it's value cannot change).
    source_set = set(variables)
    for _, change_var in petri_net.nodes(data="change"):  # type: ignore
        if change_var in source_set:
            source_set.remove(change_var)  # type: ignore[reportUnknownArgumentType] # noqa
    source_nodes: list[str] = sorted(source_set)

    ctl = _create_clingo_constraints(
        variables,
        petri_net,
        problem,
        reverse_time,
        ensure_subspace,
        avoid_subspaces,
        source_nodes,
    )

    ctl.ground()
    result = ctl.solve(yield_=True)
    if isinstance(result, SolveHandle):
        with result as iterator:
            for model in iterator:
                if not on_solution(_clingo_model_to_space(model)):
                    break
    # Else: unsat, hence we don't do anything.


def trappist(
    network: BooleanNetwork | DiGraph,
    problem: str = "min",
    reverse_time: bool = False,
    solution_limit: int | None = None,
    ensure_subspace: dict[str, int] | None = None,
    avoid_subspaces: list[dict[str, int]] | None = None,
) -> list[dict[str, int]]:
    """
    Solve the given `problem` for the given `network` using the Trappist algorithm, internally relying on the
    Python bindings of the `clingo` ASP solver.

    Arguments:
        - `network`: Can be either a `BooleanNetwork`, or a Petri net (`DiGraph`) compatible with the encoding
        in `petri_net_translation` module. The behaviour is undefined for other `DiGraph` instances.
        - `problem`: `min` minimum trap spaces; `max` maximum trap spaces; `fix` fixed points. Default: `min`.
        - `reverse_time`: If `True`, a time-reversed network should be considered. Default: `False`.
        - `solution_limit`: If given, the result is limited to the given number of solutions. Default: `None`.

    The result is a list of spaces represented as dictionaries. If you want to avoid enumerating all solutions
    explicitly as one list, you can use `trappist_async` which has a similar API but can yield solutions
    one by one.

    Finally, recall that the supplied network must have its names sanitized (see `petri_net_translation` module).
    """
    if ensure_subspace is None:
        ensure_subspace = {}
    if avoid_subspaces is None:
        avoid_subspaces = []

    results: list[dict[str, int]] = []

    def save_result(x: dict[str, int]) -> bool:
        results.append(x)
        if solution_limit is None:
            return True
        else:
            return len(results) < solution_limit

    trappist_async(
        network,
        on_solution=save_result,
        problem=problem,
        reverse_time=reverse_time,
        ensure_subspace=ensure_subspace,
        avoid_subspaces=avoid_subspaces,
    )

    return results


def _clingo_model_to_space(model: Model) -> dict[str, int]:
    space: dict[str, int] = {}
    for atom in model.symbols(atoms=True):
        atom_str = str(atom)
        (variable, is_positive) = place_to_variable(atom_str)
        # This should be prevented by the "conflic-free" property of the result,
        # but just in case.
        assert variable not in space
        # Note that this is counterintuitive but correct. If "positive" symbol
        # appears in the solution, we want to fix the value to 0. This is indeed
        # the intended behaviour of the algorithm.
        space[variable] = 0 if is_positive else 1
        space[variable] = 0 if is_positive else 1
    return space


def _create_clingo_constraints(
    variables: list[str],
    petri_net: DiGraph,
    problem: str = "min",
    reverse_time: bool = False,
    ensure_subspace: dict[str, int] | None = None,
    avoid_subspaces: list[dict[str, int]] | None = None,
    optimize_source_variables: list[str] | None = None,
) -> Control:
    """
    Translate the given Petri net (represented as a `DiGraph`; see also `petri_net_translation`
    module for details) into a logic program that solves the given problem type. This logic
    program is then added to the `Control` object provided by `clingo`.

     - The `problem` arugment specifies one of the three problem types: `min` (minimum trap spaces),
     `max` (maximum trap spaces) and `fix` (fixed-points).
     - If `reverse_time` is true, the problem is solved for a time-reversed problem.
     - Argument `ensure_subspace` is a space in which all results must be included.
     - Argument `avoid_subspaces` is a list of spaces that must be avoided by all solutions.
     - Argument `optimize_source_variables` designates variables for which a `*` solution should be
     disregarded when computing maximum trap spaces.

     Finally, note that when `ensure_subspace` or `avoid_subspaces` is included, the result is
     maximal/minimal within the resulting space of solutions, not globally. For example, if specify
     some `ensure_subspace` and `problem=max`, then the result is maximal *within* that subspace,
     not globally. Furthermore, the result can still be a *superspace* of the `avoid_subspaces`
     argument. For example, if you specify that you want to avoid a particular fixed-point,
     a globally non-minimal trap space that contains this fixed-point can be still included.
    """
    if ensure_subspace is None:
        ensure_subspace = {}
    if avoid_subspaces is None:
        avoid_subspaces = []
    if optimize_source_variables is None:
        optimize_source_variables = []

    assert problem in ["min", "max", "fix"], f"Unknown problem type: {problem}"

    dom_mod = "--dom-mod=3, 16"  # for min. trap spaces and fixed points
    if problem == "max":
        dom_mod = "--dom-mod=5, 16"  # for max. trap spaces

    # "0" specifies that all solutions should be listed (we implement the limit
    # later using callbacks).
    # TODO: Explain what remaining options mean and why we need them?
    ctl = Control(["0", "--heuristic=Domain", "--enum-mod=domRec", dom_mod])

    # Declare places and their conflicts based on network variables.
    for var_name in variables:
        p_name = variable_to_place(var_name, positive=True)
        n_name = variable_to_place(var_name, positive=False)
        # Declare a positive and negative symbol.
        ctl.add(f"{{{p_name}}}.")
        ctl.add(f"{{{n_name}}}.")
        # Assert there is no conflict (we can't have both true simultaneously).
        ctl.add(f":- {p_name}, {n_name}.")
        # If we are computing fixed points, assert also that at least one should hold.
        if problem == "fix":
            ctl.add(f"{p_name} ; {n_name}.")

    # Ensure that solutions must have desired variables fixed based on `ensure_subspace`.
    for fixed_var in ensure_subspace:
        positive = True
        if ensure_subspace[fixed_var] == 1:
            positive = False
        ctl.add(f"{variable_to_place(fixed_var, positive)}.")

    # Ensure that solutions can't have variables fixed based on either subspace in `avoid_subspaces`.
    for to_avoid in avoid_subspaces:
        fixed_list = [variable_to_place(var, (to_avoid[var] != 1)) for var in to_avoid]
        fixed_vars = ", ".join(fixed_list)
        ctl.add(f":- {fixed_vars}.")

    free_places: list[str] = []
    for node, kind in petri_net.nodes(data="kind"):  # type: ignore # noqa
        assert isinstance(node, str)
        if kind == "place":
            if place_to_variable(node)[0] not in ensure_subspace:
                free_places.append(node)
        elif kind == "transition":
            if not reverse_time:
                # Compute siphons.
                predecessors: list[str] = list(petri_net.predecessors(node))  # type: ignore # noqa
                p_disjunction = "; ".join(predecessors)
                for successor in petri_net.successors(node):  # type: ignore # noqa
                    if (
                        successor not in predecessors
                    ):  # optimize obvious tautologies # noqa
                        ctl.add(f"{p_disjunction} :- {successor}.")
            else:
                # Compute traps.
                successors = list(petri_net.successors(node))  # type: ignore # noqa
                s_disjunction = "; ".join(successors)  # type: ignore # noqa
                for predecessor in petri_net.predecessors(node):  # type: ignore # noqa
                    if predecessor not in successors:
                        ctl.add(f"{s_disjunction} :- {predecessor}.")
        else:
            raise Exception(f"Unexpected node kind: `{kind}`.")

    # For maximal trap spaces, we need an extra condition.
    if problem == "max" and len(free_places) > 0:
        # Only spaces which are not fixed but the `ensure_subspace` are considered here.
        max_condition = "; ".join(free_places)
        ctl.add(f"{max_condition}.")

        # Additionally, we require source nodes to not appear as `*` when `problem=max`, since this is
        # largely useless in practice.
        for variable in optimize_source_variables:
            if variable not in ensure_subspace:
                ctl.add(
                    f"{variable_to_place(variable, True)}; {variable_to_place(variable, False)}."
                )

    return ctl


def _clingo_model_to_fixed_point(model: Model) -> dict[str, int]:
    """
    Convert a clingo `Model` to a subspace representing a single fixed point.
    That is, the space should have all model variables fixed.

    This is a "positive" conversion, because the method for computing fixed points
    produces "positive" models (i.e. the space is represented by positive atoms
    present in the model).
    """
    space: dict[str, int] = {}

    for atom in model.symbols(atoms=True):
        atom_str = str(atom)
        (variable, is_positive) = place_to_variable(atom_str)

        # This should be prevented by the "conflic-free" property of the result,
        # but just in case.
        assert variable not in space

        # Note that this is opposite to the case of trap spaces. If "positive" symbol
        # appears in the solution, we fix the value to "1". Otherwise, we fix the value to "0".

        space[variable] = 1 if is_positive else 0

    return space


def _create_clingo_fixed_point_constraints(
    variables: list[str],
    petri_net: DiGraph,
    ensure_subspace: dict[str, int] | None = None,
    avoid_subspaces: list[dict[str, int]] | None = None,
) -> Control:
    """
    Generate the ASP characterizing all deadlocks of the Petri net (equivalently all
    fixed points of the Boolean network).
    """
    if ensure_subspace is None:
        ensure_subspace = {}
    if avoid_subspaces is None:
        avoid_subspaces = []
    dom_mod = "--dom-mod=3, 16"  # for fixed points

    # "0" specifies that all solutions should be listed (we implement the limit
    # later using callbacks).
    # TODO: Explain what remaining options mean and why we need them?
    ctl = Control(["0", "--heuristic=Domain", "--enum-mod=domRec", dom_mod])

    # Declare places and their conflicts based on network variables.
    for node in variables:
        p_name = variable_to_place(node, positive=True)
        n_name = variable_to_place(node, positive=False)

        # Declare a positive and negative symbol.
        ctl.add("base", [], f"{{{p_name}}}.")
        ctl.add("base", [], f"{{{n_name}}}.")

        # Assert there is a fixed point.
        ctl.add("base", [], f":- {p_name}, {n_name}.")
        ctl.add("base", [], f"{p_name} ; {n_name}.")

    for node, kind in petri_net.nodes(data="kind"):  # type: ignore
        if kind == "place":
            continue
        elif kind == "transition":
            preds = list(petri_net.predecessors(node))  # type: ignore

            pred_rhs = "; ".join(preds)  # type: ignore
            ctl.add("base", [], f":- {pred_rhs}.")
        else:
            raise Exception(f"Unexpected node kind: `{kind}`.")

    # Ensure that solutions must have desired variables fixed based on `ensure_subspace`.
    for fixed_var, value in ensure_subspace.items():
        place_name = variable_to_place(fixed_var, positive=bool(value))
        ctl.add("base", [], f"{place_name}.")

    # Ensure that fixed points can't lie in either subspace in `avoid_subspaces`.
    for to_avoid in avoid_subspaces:
        if len(to_avoid) > 0:
            """
            Note that this is opposite to the case of trap spaces.
            m(x) = 0 ~ place b0_x and m(x) = 1 ~ place b1_x
            """
            fixed_list = [
                variable_to_place(var, (to_avoid[var] == 1)) for var in to_avoid
            ]
            fixed_vars = ", ".join(fixed_list)
            ctl.add("base", [], f":- {fixed_vars}.")
        else:
            ctl.add("base", [], "#false.")
            break  # There is no solution and we do not need to process more.

    return ctl


def compute_fixed_point_reduced_STG_async(
    petri_net: DiGraph,
    retained_set: dict[str, int],
    on_solution: Callable[[dict[str, int]], bool],
    ensure_subspace: dict[str, int] | None = None,
    avoid_subspaces: list[dict[str, int]] | None = None,
):
    """
    The same as the `compute_fixed_point_reduced_STG`, but instead of returning a
    list of fixed-points as a result, the states are returned to the supplied
    `on_solution` callback. You can stop the enumeration by
    returning `False` from this callback.
    """
    if ensure_subspace is None:
        ensure_subspace = {}
    if avoid_subspaces is None:
        avoid_subspaces = []

    # Build a copy of the original Petri net where the
    # variables in the retained set can only change their
    # value towards the "retain value".
    reduced_petri_net: DiGraph = petri_net.copy()  # type: ignore
    for node in retained_set.keys():
        b_i = retained_set[node]
        source_place = variable_to_place(node, positive=(b_i == 1))

        preds = list(reduced_petri_net.predecessors(source_place))  # type: ignore # noqa
        succs = list(reduced_petri_net.successors(source_place))  # type: ignore # noqa

        deleted_transitions = list(set(succs) - set(preds))  # type: ignore # noqa

        for trans in deleted_transitions:  # type: ignore # noqa
            reduced_petri_net.remove_node(trans)  # type: ignore # noqa

    ctl = _create_clingo_fixed_point_constraints(
        extract_variable_names(reduced_petri_net),
        reduced_petri_net,
        ensure_subspace,
        avoid_subspaces,
    )

    ctl.ground([("base", [])])
    result = ctl.solve(yield_=True)
    if isinstance(result, SolveHandle):
        with result as iterator:
            for model in iterator:
                if not on_solution(_clingo_model_to_fixed_point(model)):
                    break
    # Else: unsat, hence we don't do anything.


def compute_fixed_point_reduced_STG(
    petri_net: DiGraph,
    retained_set: dict[str, int],
    ensure_subspace: dict[str, int] = {},
    avoid_subspaces: list[dict[str, int]] = [],
    solution_limit: int | None = None,
) -> list[dict[str, int]]:
    """
    This method computes the fixed points of the given Petri-net-encoded Boolean network.
    This makes it possible to modify the Petri net instead of re-encoding the BN repeatedly
    for multiple subsequnet queries.

    The arguments `ensure_subspace`, `avoid_subspaces`, and `solution_limit` work exactly
    the same way as in the `trappist` method. Meanwhile, the `retained_set` argument is
    applied as  a restriction on the transitions of the Petri net, forcing given variables
    to retain the specified values.
    """

    results: list[dict[str, int]] = []

    def save_result(x: dict[str, int]) -> bool:
        results.append(x)
        if solution_limit is None:
            return True
        else:
            return len(results) < solution_limit

    compute_fixed_point_reduced_STG_async(
        petri_net,
        retained_set,
        on_solution=save_result,
        ensure_subspace=ensure_subspace,
        avoid_subspaces=avoid_subspaces,
    )
    return results
