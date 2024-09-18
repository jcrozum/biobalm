from __future__ import annotations

from typing import TYPE_CHECKING, cast, Literal

if TYPE_CHECKING:
    from biobalm.succession_diagram import SuccessionDiagram
    from biodivine_aeon import VariableId

from biodivine_aeon import (
    AsynchronousGraph,
    ColoredVertexSet,
    VertexSet,
    Attractors,
    Reachability,
)
from biobalm.symbolic_utils import state_list_to_bdd
from biobalm.space_utils import is_subspace, intersect
from biobalm.types import BooleanSpace
import biodivine_aeon
import copy


def symbolic_attractor_fallback(
    sd: SuccessionDiagram,
    node_id: int,
) -> tuple[list[BooleanSpace], list[VertexSet]]:
    """
    In case the attractor candidates cannot be computed, this fallback method
    can be used to attempt fully symbolic attractor computation for the given node.

    If the normal method works, then it is usually faster, but in the rare cases
    where it fails, this can sometimes help to solve the few outlier nodes that
    would otherwise block the computation.
    """

    old_log_level = biodivine_aeon.LOG_LEVEL

    if sd.config["debug"]:
        print(f"[{node_id}] > Start symbolic fallback.")
        biodivine_aeon.LOG_LEVEL = biodivine_aeon.LOG_ESSENTIAL

    node_data = sd.node_data(node_id)
    node_space = node_data["space"]

    candidates = sd.symbolic.mk_subspace(node_space)
    if node_data["expanded"]:
        # This node has successors that we should exclude from the attractor search.
        for s in sd.node_successors(node_id, compute=False):
            s_space = sd.node_data(s)["space"]
            candidates = candidates.minus(sd.symbolic.mk_subspace(s_space))

    if node_data["skipped"]:
        # This is the same method that we applied to candidate states computation.
        initial_size = candidates.cardinality()
        for n in sd.node_ids():
            n_data = sd.node_data(n)
            if is_subspace(node_space, n_data["space"]):
                continue
            if n_data["attractor_candidates"] == [] or n_data["attractor_seeds"] == []:
                # This will create a lot of duplicates, but it seems to be better than
                # not doing it at all.
                common_subspace = intersect(node_space, n_data["space"])
                if common_subspace is not None:
                    candidates = candidates.minus(
                        sd.symbolic.mk_subspace(common_subspace)
                    )
        if sd.config["debug"]:
            print(
                f"[{node_id}] Simplified symbolic fallback candidates with skip nodes from {initial_size} to {candidates.cardinality()}."
            )

    # These should cover all cycles in the network, so transition-guided reduction
    # only needs to consider these variables.
    internal_nfvs = sd.node_percolated_nfvs(node_id, compute=True)

    if sd.config["debug"]:
        print(f"[{node_id}] > Initial attractor candidates: {candidates}")

    candidates = Attractors.transition_guided_reduction(
        sd.symbolic,
        candidates,
        internal_nfvs,
    )

    if not sd.node_is_minimal(node_id) and not candidates.is_empty():
        # This like it could be a motif-avoidant attractor. We better investigate this further,
        # because this could get complicated...
        avoid = sd.symbolic.mk_empty_colored_vertices()
        if node_data["expanded"]:
            for s in sd.node_successors(node_id):
                s_space = sd.node_data(s)["space"]
                avoid = avoid.union(sd.symbolic.mk_subspace(s_space))

        if sd.config["debug"]:
            print(
                f"[{node_id}] > Reduction was ineffective. Start computing avoid set: {avoid}"
            )

        avoid = Reachability.reach_bwd(sd.symbolic, avoid)

        if sd.config["debug"]:
            print(f"[{node_id}] > Avoid set: {avoid}")

        candidates = candidates.minus(avoid)
    else:
        if sd.config["debug"]:
            print(
                f"[{node_id}] > Node is minimal and the reduction did not finish with an empty set. This is fine."
            )

    if sd.config["debug"]:
        print(f"[{node_id}] > Attractor candidates after reduction: {candidates}")

    attractors = Attractors.xie_beerel(sd.symbolic, candidates)

    biodivine_aeon.LOG_LEVEL = old_log_level

    result_seeds: list[BooleanSpace] = []
    result_sets: list[VertexSet] = []
    for attr in attractors:
        attr_vertices = attr.vertices()
        attr_seed = next(attr_vertices.items()).to_dict()
        attr_seed_named = {
            sd.network.get_variable_name(k): v for (k, v) in attr_seed.items()
        }
        result_seeds.append(cast(BooleanSpace, attr_seed_named))
        result_sets.append(attr_vertices)

    return (result_seeds, result_sets)


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
    all_conflict_vars = copy.copy(conflict_vars)

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

    # This algorithm has two modes of operation: If `avoid` is `None`,
    # then the only thing that we can do is compute forward reachability.
    # However, if we have any avoid states, we can run two reachability
    # operations that interleave each other: One forward from the
    # pivot vertex, the other backward from the `avoid` set. If these two
    # sets ever intersect, we know that there is a path from pivot into the
    # `avoid` set. To further speed up the computation, we first prioritize
    # the variables in which the pivot differs from the avoid set (so called
    # conflict variables). If all such variables are covered, we then consider
    # those that are close to the conflict variables in terms of regulations
    # (i.e. if we can't update a conflict variable, we ideally want to update
    # one that regulates it).

    # True if the main cycle completes with all saturation procedures fully
    # completed and no unprocessed variables remaining.
    all_done = False

    while not all_done:
        all_done = True

        # Saturate reach_set with currently selected variables, but only if
        # it's symbolic size is smaller than that of the avoid set (reach set
        # tends to grow quite large and we'd like to avoid that).
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
                    all_done = False  # The main loop should continue.
                    updated = reach_set.union(successors)
                    no_avoid = avoid is None
                    avoid_is_larger = (
                        avoid is not None
                    ) and avoid.symbolic_size() >= updated.symbolic_size()
                    all_variables_done = (
                        len(conflict_vars) == 0 and len(other_vars) == 0
                    )
                    if no_avoid or avoid_is_larger or all_variables_done:
                        reach_set = updated
                        saturation_done = False
                        if reach_set.symbolic_size() > 100_000 and sd.config["debug"]:
                            print(
                                f"[{node_id}] > Saturation({len(saturated_vars)}) Incremented forward reach set: {reach_set}"
                            )
                        break

        if avoid is not None:
            # If `avoid` is not `None`, we also want to expand it backwards.
            saturation_done = False
            while not saturation_done:
                if not avoid.intersect(reach_set).is_empty():
                    if sd.config["debug"]:
                        print(f"[{node_id}] > Discovered avoid state. Done.")
                    return None

                saturation_done = True
                for var in saturated_vars:
                    predecessors = graph.var_pre_out(var, avoid)
                    if not predecessors.is_empty():
                        all_done = False
                        avoid = avoid.union(predecessors)
                        saturation_done = False
                        if avoid.symbolic_size() > 100_000 and sd.config["debug"]:
                            print(
                                f"[{node_id}] > Saturation({len(saturated_vars)}) Incremented backward avoid set: {avoid}"
                            )
                        break

        if sd.config["debug"] and (
            reach_set.symbolic_size() > 100_000
            or (avoid is not None and avoid.symbolic_size() > 100_000)
        ):
            print(
                f"[{node_id}] > Saturation({len(saturated_vars)}) Finished with avoid set {avoid} and reach set {reach_set}."
            )

        # Once saturated, try to expand the saturated
        # collection with either a conflict variable or
        # other variable.

        # First, we sort the non-conflict variables based on their
        # distance towards the conflict variables.
        network = sd.node_percolated_network(node_id)
        distances: dict[VariableId, int] = {
            v: network.variable_count() for v in network.variables()
        }
        visited: set[VariableId]
        if len(conflict_vars) > 0:
            visited = set(conflict_vars)
            current_level = set(conflict_vars)
        else:
            visited = set(all_conflict_vars)
            current_level = set(all_conflict_vars)

        next_level: set[VariableId] = set()

        distance = 0
        while len(current_level) > 0:
            for var in current_level:
                distances[var] = min(distance, distances[var])
                for s in network.predecessors(var):
                    if s not in visited:
                        visited.add(s)
                        next_level.add(s)
            current_level = next_level
            next_level = set()
            distance += 1

        other_sorted = sorted(other_vars, key=lambda x: distances[x])
        for var in conflict_vars + other_sorted:
            can_go_fwd = graph.var_post_out(var, reach_set)
            if avoid is not None:
                can_go_bwd = graph.var_pre_out(var, avoid)
            else:
                can_go_bwd = graph.mk_empty_colored_vertices()

            if can_go_fwd.is_empty() and can_go_bwd.is_empty():
                continue

            all_done = False

            reach_set = reach_set.union(can_go_fwd)
            if avoid is not None:
                avoid = avoid.union(can_go_bwd)

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
