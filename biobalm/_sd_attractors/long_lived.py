from typing import cast
from copy import copy
from biobalm.succession_diagram import SuccessionDiagram
from biobalm.types import BooleanSpace
from biodivine_aeon import AsynchronousGraph, ColoredVertexSet


def compute_long_lived_phenotypes(
    sd: SuccessionDiagram,
    node_id: int,
    only_maximal: bool = True,
) -> list[BooleanSpace]:
    """
    Compute the phenotypes defined by long-lived components present
    in the given succession diagram node.

    Parameters
    ----------
    sd : SuccessionDiagram
        The SuccessionDiagram instance containing the state space.
    node_id : int
        The ID of the node for which to compute long-lived phenotypes.

    Returns
    -------
    list[BooleanSpace]
        A list of BooleanSpace objects representing the long-lived phenotypes.
    """
    node_data = sd.node_data(node_id)

    if not node_data["expanded"]:
        # TODO: We could compute long-lived phenotypes without expanding the node,
        # but this would include extra phenotypes from the smaller trap spaces.
        raise RuntimeError(
            f"Node {node_id} is not expanded. Cannot compute long-lived phenotypes yet."
        )

    node_network = sd.node_percolated_network(node_id, compute=True)

    if node_network.variable_count() == 0:
        # This is a fixed-point.
        return [copy(node_data["space"])]

    stg = AsynchronousGraph(node_network)

    # Just checking that the network does not have any parameters at this point.
    assert stg.mk_unit_colors().is_singleton()

    search_space = stg.mk_unit_colored_vertices()

    for successor in sd.node_successors(node_id):
        motif = sd.edge_stable_motif(node_id, successor, reduced=True)
        motif_space = stg.mk_subspace(motif)
        search_space = search_space.minus(motif_space)

    phenotypes: set[str] = set()

    to_explore: list[tuple[ColoredVertexSet, BooleanSpace]] = [(search_space, {})]

    while len(to_explore) > 0:
        min_space_size = min([len(x[1]) for x in to_explore])
        to_search = [x for x in to_explore if len(x[1]) == min_space_size]
        to_explore = [x for x in to_explore if len(x[1]) > min_space_size]

        print(
            f"Searching for maximal components in {len(to_search)} spaces (fixed vars: {min_space_size}/{node_network.variable_count()})."
        )

        long_lived_components: set[ColoredVertexSet] = set()

        for i, (space, phenotype) in enumerate(to_search):
            print(
                f"Searching space {i+1}/{len(to_search)} with {min_space_size} fixed vars."
            )
            for component in maximal_long_lived_components(stg, space):
                if component not in long_lived_components:
                    long_lived_components.add(component)
                    print(
                        f"Found a long-lived component ({len(long_lived_components)} found so far)"
                    )
                # assert component.is_subset(space)
                # component_subspace = component.vertices().enclosing_named_subspace()
                # assert component_subspace is not None
                # assert len(component_subspace) >= len(phenotype), f"{component_subspace} >= {phenotype} {subspace(stg, component)} {subspace(stg, space)}"

        for i, component in enumerate(long_lived_components):
            inferred_phenotype = component.vertices().enclosing_named_subspace()
            assert inferred_phenotype is not None
            phenotype = {k: 1 if v else 0 for k, v in inferred_phenotype.items()}
            assert phenotype is not None
            phenotype_str = str(sorted(phenotype.items()))
            if phenotype_str not in phenotypes:
                phenotypes.add(phenotype_str)
                print(
                    f"[{i+1} / {len(long_lived_components)}] Found a new phenotype with {len(phenotype)} fixed vars ({len(phenotypes)} found so far)."
                )
            else:
                print(
                    f"[{i+1} / {len(long_lived_components)}] Duplicate phenotype (size {len(phenotype)})."
                )

            if not only_maximal:
                # Enqueue all smaller subspace combinations
                for var in stg.network_variable_names():
                    if var in phenotype:
                        continue  # This variable is fixed by the SCC.
                    var_true: ColoredVertexSet = trim(
                        stg, component.intersect(stg.mk_subspace({var: True}))
                    )
                    var_false: ColoredVertexSet = trim(
                        stg, component.intersect(stg.mk_subspace({var: False}))
                    )
                    if not var_true.is_empty() and not has_percolation(stg, var_true):
                        true_subspace = var_true.vertices().enclosing_named_subspace()
                        assert true_subspace is not None
                        found = False
                        # We have to check this explicitly, because it might have gotten there also from some larger value
                        # of min_space_size, not just the one we are computing now.
                        for x, _ in to_explore:
                            if x.is_subset(var_true) and var_true.is_subset(x):
                                found = True
                                break
                        if not found:
                            print(
                                f"\t\tNew search space with {len(true_subspace)} fixed vars ({var} = 1)."
                            )
                            to_explore.append(
                                (
                                    var_true,
                                    {
                                        k: 1 if v else 0
                                        for k, v in true_subspace.items()
                                    },
                                )
                            )
                    if not var_false.is_empty() and not has_percolation(stg, var_false):
                        false_subspace = var_false.vertices().enclosing_named_subspace()
                        assert false_subspace is not None
                        print(
                            f"\t\tNew search space with {len(false_subspace)} fixed vars ({var} = 0)."
                        )
                        found = False
                        for x, _ in to_explore:
                            if x.is_subset(var_false) and var_false.is_subset(x):
                                found = True
                                break
                        if not found:
                            to_explore.append(
                                (
                                    var_false,
                                    {
                                        k: 1 if v else 0
                                        for k, v in false_subspace.items()
                                    },
                                )
                            )

    return cast(
        list[BooleanSpace],
        [(dict(eval(x)) | node_data["space"]) for x in sorted(phenotypes)],
    )


def subspace(stg: AsynchronousGraph, data: ColoredVertexSet) -> BooleanSpace:
    result: BooleanSpace = {}
    for v in stg.network_variable_names():
        v_true = data.intersect(stg.mk_subspace({v: True}))
        v_false = data.intersect(stg.mk_subspace({v: False}))
        if v_true.is_empty():
            result[v] = 0
        elif v_false.is_empty():
            result[v] = 1
    return result


def maximal_long_lived_components(
    stg: AsynchronousGraph, search_space: ColoredVertexSet
) -> list[ColoredVertexSet]:
    search_space = trim(stg, search_space)
    if search_space.is_empty() or has_percolation(stg, search_space):
        return []

    # Invariant: Every set inserted to this queue must:
    #  - be non-empty
    #  - be non-percolable
    #  - be trimmed
    to_explore = [search_space]

    result: list[ColoredVertexSet] = []

    while len(to_explore) > 0:
        universe = to_explore.pop()

        # print(f"Start processing universe: {universe}. Remaining: {sum([ x.cardinality() for x in to_explore ])} states in {len(to_explore)} sets.")

        pivot = universe.pick_vertex()

        bwd = reach_bwd(stg, pivot, universe)
        scc = reach_fwd(stg, pivot, bwd)

        bwd_rest = trim(stg, bwd.minus(scc))
        universe_rest = trim(stg, universe.minus(bwd))

        if not bwd_rest.is_empty() and not has_percolation(stg, bwd_rest):
            to_explore.append(bwd_rest)
        if not universe_rest.is_empty() and not has_percolation(stg, universe_rest):
            to_explore.append(universe_rest)

        if not has_percolation(stg, scc):
            print(f"Found long-lived SCC {scc}.")
            result.append(scc)
        # else:
        #    print(f"Skipping short-lived SCC {scc}.")

    return result


def has_percolation(stg: AsynchronousGraph, data: ColoredVertexSet) -> bool:
    for v in stg.network_variables():
        can_escape = stg.var_can_post_out(v, data)
        if data.is_subset(can_escape):
            return True
    return False


def trim(stg: AsynchronousGraph, data: ColoredVertexSet) -> ColoredVertexSet:
    # print(f"Start trimming set {data}.")
    while True:
        has_successor = stg.pre(data).intersect(data)
        has_predecessor = stg.post(data).intersect(data)
        non_trim = has_successor.intersect(has_predecessor)
        # print(f"Trimming to set {non_trim}.")
        if data.is_subset(non_trim):
            return data
        data = non_trim


def reach_fwd(
    stg: AsynchronousGraph, initial: ColoredVertexSet, universe: ColoredVertexSet
) -> ColoredVertexSet:
    # print(f"Start forward reachability for set {initial} in {universe}.")
    result = initial
    while True:
        done = True
        for var in reversed(stg.network_variables()):
            next = stg.var_post(var, result).intersect(universe)
            if not next.is_subset(result):
                result = result.union(next)
                # print(f"Reachability increased to {result} via {var}.")
                done = False
                break
        if done:
            return result


def reach_bwd(
    stg: AsynchronousGraph, initial: ColoredVertexSet, universe: ColoredVertexSet
) -> ColoredVertexSet:
    # print(f"Start backward reachability for set {initial} in {universe}.")
    result = initial
    while True:
        done = True
        for var in reversed(stg.network_variables()):
            prev = stg.var_pre(var, result).intersect(universe)
            if not prev.is_subset(result):
                result = result.union(prev)
                # print(f"Reachability increased to {result} via {var}.")
                done = False
                break
        if done:
            return result
