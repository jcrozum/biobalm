from __future__ import annotations

import itertools as it
import sys
from typing import TYPE_CHECKING, Callable, cast

from biodivine_aeon import BooleanNetwork, SymbolicContext

from biobalm._sd_algorithms.expand_bfs import expand_bfs
from biobalm.space_utils import percolate_network, percolate_space
from biobalm.types import BooleanSpace

if TYPE_CHECKING:
    from biobalm.succession_diagram import SuccessionDiagram

    ExpanderFunctionType = Callable[
        [SuccessionDiagram, int | None, int | None, int | None],
        bool,
    ]

sys.path.append(".")

DEBUG = False


def expand_source_SCCs(
    sd: SuccessionDiagram,
    expander: ExpanderFunctionType = expand_bfs,
    check_maa: bool = True,
) -> bool:
    """
    1. percolate
    2. find source nodes, fix combinations and then percolate
    3. at every level, find all source SCCs, expand them
    4. when there are no more SCCs, expand in a usual way.

    TODO: the function uses bfs for source SCCs and for remaining nodes.
    it should be able to import any expansion strategy and use it.

    TODO: it will be rare especially for empirical models,
    but there can be source SCCs even before fixing source node combinations.
    It could be useful to find them first, calculate scc_sd and store them somewhere,
    rather than calculating the same thing for every source node combination

    Parameters
    ----------
    check_maa - if True, check for motif avoidant attractors when
                expanding the sd with source SCCs.
                This grants significant speedup when searching for all attractors.
                if False, ignores motif avoidant attractors in source SCCs.

    """

    root = sd.root()

    current_level = [root]
    next_level: list[int] = []
    final_level: list[int] = []  # from here there are no more source SCCs

    # percolate constant nodes
    perc_space = percolate_space(sd.symbolic, {})
    sd.node_data(root)["space"] = perc_space

    # find source nodes
    perc_bn = percolate_network(sd.network, perc_space)
    source_nodes = find_source_nodes(perc_bn)

    # get source nodes combinations and expand root node
    if len(source_nodes) != 0:
        bin_values_iter = it.product(range(2), repeat=len(source_nodes))
        for bin_values in bin_values_iter:
            source_comb = cast(BooleanSpace, dict(zip(source_nodes, bin_values)))

            sub_space = source_comb
            sub_space.update(perc_space)

            next_level.append(sd._ensure_node(root, sub_space))  # type: ignore

        sd.node_data(root)["expanded"] = True
        sd.node_data(root)[
            "attractor_seeds"
        ] = []  # no need to look for attractors here
        current_level = next_level
        next_level = []

    while len(current_level) > 0:
        if DEBUG:
            print(f"{current_level=}")

        # each level consists of one round of fixing all source SCCs
        for node_id in current_level:
            sub_sds = list(sd.component_subdiagrams(node_id))
            # if there are no more source SCCs in this node, move it to the final level
            if not sub_sds:
                final_level.append(node_id)
                continue

            # attach all source scc sd one by one
            current_branches = [
                node_id
            ]  # this is where the scc_sd should be "attached"
            next_branches: list[int] = []
            for sub_network_diagram in sub_sds:
                scc_sd, exist_maa = find_subnetwork_sd(
                    sub_network_diagram, expander, check_maa
                )

                if exist_maa:  # we check for maa, and it exists
                    continue
                # add scc_sd to every single branch in the original sd
                for branch in current_branches:
                    new_branches = attach_scc_sd(sd, scc_sd, branch, check_maa)
                    next_branches.extend(new_branches)

                current_branches = next_branches
                next_branches = []

            # nothing was attached
            # This happens when there were only source SCCs with motif avoidant attractors
            if current_branches == [node_id]:
                final_level.append(node_id)
            # once all the source SCCs are expanded, move on to the next level
            # the final branches become the next level
            else:
                next_level.extend(current_branches)

        current_level = next_level
        next_level = []

    # Finally, expand the sd in a usual way.
    # TODO: motif avoidance that is already calculated in the source SCCs should be used here somehow.
    if DEBUG:
        print(f"{final_level=}")
    for node_id in final_level:
        # These assertions should be unnecessary, but just to be sure.
        assert not sd.node_data(node_id)["expanded"]  # expand nodes from here
        assert (
            sd.node_data(node_id)["attractor_seeds"] is None
        )  # check attractors from here

        # restore this once we allow all expansion algorithms to expand from a node
        # expander(sd, node_id)
        sd.expand_bfs(node_id)

    return True


def find_source_nodes(
    network: BooleanNetwork, ctx: SymbolicContext | None = None
) -> list[str]:
    """
    Return the "source nodes" of a `BooleanNetwork`. That is, variables whose value
    cannot change, but is not fixed to a `true`/`false` constant.

    Note that this internally uses BDD translation to detect identity functions
    semantically rather than syntactically. If you already have a `SymbolicContext`
    for the given `network` available, you can supply it as the second argument.
    """
    if ctx is None:
        ctx = SymbolicContext(network)

    result: list[str] = []
    for var in network.variables():
        update_function = network.get_update_function(var)
        assert update_function is not None  # Parameters are not allowed.
        fn_bdd = ctx.mk_update_function(update_function)
        var_bdd = ctx.mk_network_variable(var)
        if fn_bdd == var_bdd:
            result.append(network.get_variable_name(var))

    return result


def find_subnetwork_sd(
    sub_network_diagram: SuccessionDiagram,
    expander: ExpanderFunctionType,
    check_maa: bool,
) -> tuple[SuccessionDiagram, bool]:
    """
    Computes a `SuccessionDiagram` of a particular sub-network using an expander function.

    If `check_maa` is set to `True`, also checks if the resulting succession diagram admits
    motif avoidant attractors.

    Returns
    -------
    scc_sd - succession diagram of the source SCC
    exist_maa - False if motif avoidance is not checked, or if motif avoidance does not exist
                True if there is motif avoidance

    """
    fully_expanded = expander(sub_network_diagram, None, None, None)
    assert fully_expanded

    has_maa = False
    if check_maa:
        # check for motif avoidant attractors
        # TODO: somehow skip this calculation when this source SCC appears again later.
        # it will appear again, since souce SCC with maa are not fixed.
        motif_avoidant_count = 0
        for node in sub_network_diagram.node_ids():
            attr = sub_network_diagram.node_attractor_seeds(node, compute=True)
            if not sub_network_diagram.node_is_minimal(node):
                motif_avoidant_count += len(attr)
        if motif_avoidant_count != 0:
            # ignore source SCCs with motif avoidant attractors
            has_maa = True

    return sub_network_diagram, has_maa


def attach_scc_sd(
    sd: SuccessionDiagram,
    scc_sd: SuccessionDiagram,
    branch: int,
    check_maa: bool,
) -> list[int]:
    """
    Attach scc_sd to the given branch point of the sd.
    Returns the new branching points.

    Parameters
    ----------
    check_maa - If true, it is assumed that scc_sd is already checked
                to not have motif avoidant attractors.
                Make sure only give scc_sd with no maa
    """
    if len(scc_sd) == 1:
        return [branch]

    next_branches: list[int] = []
    size_before_attach = sd.dag.number_of_nodes()
    # first add all the nodes using their first parent
    for scc_node_id in scc_sd.node_ids():
        if scc_node_id == 0:  # no need to add the root of scc_sd
            continue

        scc_parent_id = cast(
            int,
            list(scc_sd.dag.predecessors(scc_node_id))[0],  # type: ignore
        )  # get the first parent
        assert scc_parent_id < scc_node_id

        if scc_parent_id == 0:  # the parent is the root of scc_sd
            parent_id = branch
        else:
            parent_id = size_before_attach + scc_parent_id - 1

        motif = scc_sd.edge_stable_motif(scc_parent_id, scc_node_id)
        motif.update(cast(BooleanSpace, sd.dag.nodes[branch]["space"]))

        child_id = sd._ensure_node(parent_id, motif)  # type: ignore
        if check_maa:
            sd.dag.nodes[parent_id][
                "attractors"
            ] = []  # no need to check for attractors in these nodes

        if scc_sd.node_is_minimal(scc_node_id) and child_id not in next_branches:
            next_branches.append(child_id)

    # now add all the missing edges
    for scc_node_id in scc_sd.node_ids():
        if scc_node_id == 0:
            parent_id = branch
        else:
            parent_id = size_before_attach + scc_node_id - 1

        scc_child_ids = cast(list[int], list(scc_sd.dag.successors(scc_node_id)))  # type: ignore
        for scc_child_id in scc_child_ids:
            motif = scc_sd.edge_stable_motif(scc_node_id, scc_child_id)
            motif.update(cast(BooleanSpace, sd.dag.nodes[branch]["space"]))

            child_id = sd._ensure_node(parent_id, motif)  # type: ignore
            assert child_id == size_before_attach + scc_child_id - 1

        # if the node had any child node, consider it expanded.
        if len(scc_child_ids) > 0:
            sd.dag.nodes[parent_id]["expanded"] = True

    return next_branches
