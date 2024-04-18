from __future__ import annotations

from typing import TYPE_CHECKING, cast, Literal

if TYPE_CHECKING:
    from biobalm.succession_diagram import SuccessionDiagram
    from biodivine_aeon import VariableId

from biodivine_aeon import AsynchronousGraph, ColoredVertexSet, VertexSet
from biobalm.symbolic_utils import state_list_to_bdd
from biobalm.types import BooleanSpace


def compute_attractors_symbolic(
    sd: SuccessionDiagram,
    node_id: int,
    candidate_states: list[BooleanSpace],
    seeds_only: bool = False,
) -> tuple[list[BooleanSpace], list[VertexSet] | None]:
    """
    Uses exhaustive symbolic reachability to eliminate spurious candidate states
    and to compute the exact attractor sets.

    Parameters
    ----------
    sd : SuccessionDiagram
        The succession diagram in question.
    node_id : int
        An ID of an SD node for which we are computing the candidates.
    candidate_states: bool
        The list of candidate states covering the attractors in the
        given node.
    seeds_only: bool
        If `True`, the method can terminate early once it is guaranteed that
        all seeds have been identified. In such case, the list of returned sets
        is `None`.

    Returns
    -------
    tuple[list[BooleanSpace], list[VertexSet] | None]
        The list of attractor seed states and corresponding attractor sets (if computed).
    """

    node_data = sd.node_data(node_id)
    node_space = node_data["space"]
    bn_reduced = sd.node_percolated_network(node_id, compute=True)
    graph_reduced = AsynchronousGraph(bn_reduced)
    symbolic_ctx = graph_reduced.symbolic_context()

    # Ignore variables that are already fixed in the nodes space.
    # These are not used in the reduced network.
    candidate_states_reduced: list[BooleanSpace] = []
    for candidate in candidate_states:
        candidate_reduced: BooleanSpace = {
            k: v for (k, v) in candidate.items() if k not in node_space
        }
        candidate_states_reduced.append(candidate_reduced)

    child_motifs_reduced = []
    if node_data["expanded"]:
        children = sd.node_successors(node_id, compute=False)
        child_motifs_reduced = [
            sd.edge_stable_motif(node_id, s, reduced=True) for s in children
        ]

    child_motifs_bdd = state_list_to_bdd(symbolic_ctx, child_motifs_reduced)
    candidate_bdd = state_list_to_bdd(symbolic_ctx, candidate_states_reduced)

    children_set = ColoredVertexSet(symbolic_ctx, child_motifs_bdd)
    candidate_set = ColoredVertexSet(symbolic_ctx, candidate_bdd)

    avoid = candidate_set.union(children_set)

    if sd.config["debug"]:
        print(
            f"[{node_id}] > Start symbolic seed state identification with {len(candidate_states)} candidates and {avoid} avoid states."
        )

    seeds: list[BooleanSpace] = []
    sets: list[ColoredVertexSet] = []
    for i, candidate in enumerate(candidate_states_reduced):
        is_last = i == len(candidate_states_reduced) - 1
        is_minimal = len(child_motifs_reduced) == 0
        if seeds_only and is_minimal and is_last and len(seeds) == 0:
            # This node is pseudo-minimal, so it must contain at least
            # one attractor seed. If we only care about the attractor seeds,
            # we can thus return the last candidate without checking it.
            if sd.config["debug"]:
                print(
                    f"[{node_id}] > Single seed remaining in a (pseduo) minimal space. Done."
                )
            return ([candidate | node_space], None)
        candidate_singleton = graph_reduced.mk_subspace(candidate)

        avoid = avoid.minus(candidate_singleton)

        closure = symbolic_attractor_test(sd, node_id, graph_reduced, candidate, avoid)

        if closure is None:
            # This candidate can reach someone else in the candidate set,
            # or they can reach one of the child spaces. Hence it is not
            # an attractor.
            continue

        # Otherwise, we have an attractor set and we can save it.
        # It also becomes a part of the avoid set.
        avoid = avoid.union(closure)
        seeds.append(candidate | node_space)
        sets.append(closure)

    if sd.config["debug"]:
        print(f"[{node_id}] > Finished identification with {len(seeds)} seed states.")

    space_symbolic = sd.symbolic.mk_subspace(node_space).vertices()
    sets_converted: list[VertexSet] = []
    for s in sets:
        # Extend the attractor set with fixed nodes from the node space.
        vertices = s.vertices()
        vertices = sd.symbolic.transfer_from(vertices, graph_reduced)
        vertices = vertices.intersect(space_symbolic)
        sets_converted.append(vertices)

    return (seeds, sets_converted)


def symbolic_attractor_test(
    sd: SuccessionDiagram,
    node_id: int,
    graph: AsynchronousGraph,
    pivot: BooleanSpace,
    avoid_set: ColoredVertexSet,
) -> ColoredVertexSet | None:
    """
    Use symbolic reachability with saturation to compute the set of states
    reachable from the given `pivot`, or `None` if the reachability procedure
    intersects with any state in the `avoid` set.

    *The reason why we use `ColoredVertexSet` instead of `VertexSet` is mostly
    a technicality that should be irelevant in biobalm, since we don't allow any
    parameters outside of unknown inputs.*
    """

    if avoid_set.is_empty():
        avoid = None
    else:
        avoid = avoid_set

    reach_set = graph.mk_subspace(pivot)

    # The variables for which we already maintain that `reach_set`
    # contains all reachable states.
    saturated_vars: list[VariableId] = []

    # Variables where `reach_set` differs from the states in `avoid`.
    # We should prioritize updating these variables, because they *need*
    # to be updated if we are ever to reach `avoid`.
    conflict_vars: list[VariableId] = []

    # Populate conflict vars, assuming we have any.
    if avoid is not None:
        avoid_bdd = avoid.to_bdd()
        for name, val in pivot.items():
            var = graph.find_network_variable(name)
            bdd_var = graph.symbolic_context().find_network_bdd_variable(name)
            assert var is not None and bdd_var is not None
            incompatible = avoid_bdd.r_select({bdd_var: cast(Literal[0, 1], 1 - val)})
            if not incompatible.is_false():
                conflict_vars.append(var)
        conflict_vars = sort_variable_list(conflict_vars)

    # Remaining network variables that are still relevant, but may not
    # be necessary to reach `avoid`.
    other_vars: list[VariableId] = [
        x for x in graph.network_variables() if x not in conflict_vars
    ]
    other_vars = sort_variable_list(other_vars)

    if sd.config["debug"]:
        print(
            f"[{node_id}] > Start symbolic reachability with {len(conflict_vars)} conflict variables and {len(other_vars)} other variables."
        )

    all_done = False
    while not all_done:
        all_done = True

        # Saturate reach set with currently selected variables.
        saturation_done = False
        while not saturation_done:
            if avoid is not None and not avoid.intersect(reach_set).is_empty():
                if sd.config["debug"]:
                    print(f"[{node_id}] > Discovered avoid state. Done.")
                return None

            saturation_done = True
            for var in saturated_vars:
                successors = graph.var_post_out(var, reach_set)
                if not successors.is_empty():
                    reach_set = reach_set.union(successors)
                    saturation_done = False
                    if reach_set.symbolic_size() > 100_000 and sd.config["debug"]:
                        print(
                            f"[{node_id}] > Saturation({len(saturated_vars)}) Expanded reach_set: {reach_set}"
                        )
                    break

        # Once saturated, try to expand the saturated
        # collection with either a conflict variable or
        # other variable.

        # First try conflict vars, then other vars.
        for var in conflict_vars + other_vars:
            successors = graph.var_post_out(var, reach_set)
            if not successors.is_empty():
                reach_set = reach_set.union(successors)
                all_done = False

                # This is a bit wasteful but at this point it
                # should be irrelevant for performance.
                if var in conflict_vars:
                    conflict_vars.remove(var)

                if var in other_vars:
                    other_vars.remove(var)

                saturated_vars.append(var)
                saturated_vars = sort_variable_list(saturated_vars)

                if sd.config["debug"]:
                    print(
                        f"[{node_id}] > Saturation({len(saturated_vars)}) Added saturation variable. {len(conflict_vars)} conflict and {len(other_vars)} other variables remaining."
                    )
                break

    if sd.config["debug"]:
        print(f"[{node_id}] > Reachability completed with {reach_set}.")

    return reach_set


def sort_variable_list(variables: list[VariableId]):
    return list(sorted(variables, reverse=True))
