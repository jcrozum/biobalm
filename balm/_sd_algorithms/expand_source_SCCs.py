from __future__ import annotations

import itertools as it
import sys
from typing import TYPE_CHECKING, Callable, Iterable, cast

import networkx as nx  # type: ignore
from biodivine_aeon import BooleanNetwork
from networkx import DiGraph

import balm.SuccessionDiagram
from balm._sd_algorithms.expand_bfs import expand_bfs
from balm.interaction_graph_utils import infer_signed_interaction_graph
from balm.petri_net_translation import extract_variable_names, network_to_petrinet
from balm.space_utils import percolate_network, percolate_space

if TYPE_CHECKING:
    expander_function_type = Callable[
        [balm.SuccessionDiagram.SuccessionDiagram, int | None, int | None, int | None],
        bool,
    ]

sys.path.append(".")

DEBUG = False


def expand_source_SCCs(
    sd: balm.SuccessionDiagram.SuccessionDiagram,
    expander: expander_function_type = expand_bfs,
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
    perc_space = percolate_space(sd.network, {}, strict_percolation=False)
    sd.dag.nodes[root]["space"] = perc_space

    # find source nodes
    perc_bn = percolate_network(sd.network, perc_space)
    source_nodes = find_source_nodes(perc_bn)

    # get source nodes combinations and expand root node
    if len(source_nodes) != 0:
        bin_values_iter = it.product(range(2), repeat=len(source_nodes))
        for bin_values in bin_values_iter:
            source_comb = dict(zip(source_nodes, bin_values))

            sub_space = source_comb
            sub_space.update(perc_space)

            next_level.append(sd._ensure_node(root, sub_space))  # type: ignore

        sd.dag.nodes[root]["expanded"] = True
        sd.dag.nodes[root]["attractors"] = []  # no need to look for attractors here
        current_level = next_level
        next_level = []

    while len(current_level) > 0:
        if DEBUG:
            print(f"{current_level=}")

        # each level consists of one round of fixing all source SCCs
        for node_id in current_level:
            sub_space = cast(dict[str, int], sd.dag.nodes[node_id]["space"])

            # find source SCCs
            clean_bnet, clean_bn = perc_and_remove_constants_from_bn(perc_bn, sub_space)
            source_scc_list = find_source_SCCs(clean_bn)
            if DEBUG:
                print(f"{source_scc_list=}")

            # if there are no more source SCCs in this node, move it to the final level
            if len(source_scc_list) == 0:
                final_level.append(node_id)
                continue

            # attach all source scc sd one by one
            current_branches = [
                node_id
            ]  # this is where the scc_sd should be "attached"
            next_branches: list[int] = []
            while len(source_scc_list) > 0:
                source_scc = source_scc_list.pop(0)
                scc_sd, exist_maa = find_scc_sd(
                    clean_bnet, source_scc, expander, check_maa
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
        assert not sd.dag.nodes[node_id]["expanded"]  # expand nodes from here
        assert sd.dag.nodes[node_id]["attractors"] is None  # check attractors from here

        # restore this once we allow all expansion algorithms to expand from a node
        # expander(sd, node_id)
        sd.expand_bfs(node_id)

    return True


def find_source_nodes(network: BooleanNetwork | DiGraph) -> list[str]:
    """
    TODO: current method translates percolated bn into petri_net.
    This is a repetition since bn is translated into petri_net when initiating an sd.
    But I don't know how to get a percolated petri_net from the original one.
    """
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

    return source_nodes


def perc_and_remove_constants_from_bn(
    bn: BooleanNetwork, space: dict[str, int]
) -> tuple[str, BooleanNetwork]:
    """
    Take a BooleanNetwork and percolate given space.
    Then remove constant nodes with rules
    A* = (A | !A)
    B* = (B & !B)
    and return a clean BooleanNetwork with no constant nodes.

    TODO: if possible, delete constant nodes from the BooleanNetwork directly,
    without converting into bnet and back.
    For now getting the percolated bnet is useful, as it is used for find_scc_sd()
    """

    perc_space = percolate_space(bn, space, strict_percolation=False)
    perc_bn = percolate_network(bn, perc_space)

    perc_bnet = perc_bn.to_bnet()

    # remove constant nodes from the bnet
    clean_bnet = ""
    for line in perc_bnet.split("\n"):
        if line.strip() == "":
            continue
        node = line.split(",")[0]
        function = line.split(",")[1].strip()
        if function == "(" + node + " | !" + node + ")":
            continue
        elif function == "(" + node + " & !" + node + ")":
            continue
        else:
            clean_bnet += line + "\n"

    clean_bn = BooleanNetwork.from_bnet(clean_bnet)

    return clean_bnet, clean_bn


def find_source_SCCs(bn: BooleanNetwork) -> list[list[str]]:
    """
    Find source SCCs given a bn.
    Note that the provided bn should not have constant nodes such as
    A* = (A | !A)
    B* = (B & !B)
    as they will show up as source SCCs.

    TODO: better way of finding source SCC,
    that does not involve translating into networkx digraph would be useful
    """
    di = infer_signed_interaction_graph(bn)
    scc_list: list[list[str]] = [
        sorted(x) for x in nx.strongly_connected_components(di)  # type: ignore
    ]

    # find source SCCs
    source_scc_list: list[list[str]] = []
    for scc in scc_list:
        source = True
        for node in scc:
            for reg in cast(Iterable[str], di.predecessors(node)):  # type: ignore
                if reg not in scc:
                    source = False
                    break
            if not source:
                break
        if source:
            source_scc_list.append(scc)

    return sorted(source_scc_list)


def find_scc_sd(
    bnet: str, source_scc: list[str], expander: expander_function_type, check_maa: bool
) -> tuple[balm.SuccessionDiagram.SuccessionDiagram, bool]:
    """
    TODO: better way that does not use bnet but rather bn directly or petri_net directly
    to find the scc_sd would be useful.
    TODO: Somehow remove implicit parameters elegantly.
    Right now, the implicit parameters are deleted from the subspaces of the returned scc_sd,
    but not from their edges.

    Returns
    -------
    scc_sd - succession diagram of the source SCC
    exist_maa - False if motif avoidance is not checked, or if motif avoidance does not exist
                True if there is motif avoidance

    """

    # get rules for only the source SCC
    scc_bnet = "targets,factors\n"
    for line in bnet.split("\n"):
        for node_name in source_scc:
            if line.startswith(node_name + ","):
                scc_bnet += line + "\n"

    scc_bn = BooleanNetwork.from_bnet(scc_bnet)

    if DEBUG:
        print("scc_bnet\n", scc_bn.to_bnet())

    # scc_bn = scc_bn.infer_regulatory_graph()

    # quick way to ignore implicit parameters coming from fake edges into source SCCs
    # e.g. in a rule such as A* = A | A & B, B may appear in the scc_bnet implicitly but is meaningless
    implicit_parameters: list[str] = []
    for var in scc_bn.implicit_parameters():
        scc_bn.set_update_function(var, "true")
        implicit_parameters.append(scc_bn.get_variable_name(var))
    scc_bn = scc_bn.infer_regulatory_graph()

    # Compute the succession diagram.
    scc_sd = balm.SuccessionDiagram.SuccessionDiagram(scc_bn)
    fully_expanded = expander(scc_sd)  # type: ignore
    assert fully_expanded

    exist_maa = False
    if check_maa:
        # check for motif avoidant attractors
        # TODO: somehow skip this calculation when this source SCC appears again later.
        # it will appear again, since souce SCC with maa are not fixed.
        motif_avoidant_count = 0
        for node in scc_sd.node_ids():
            attr = scc_sd.node_attractor_seeds(node, compute=True)
            if not scc_sd.node_is_minimal(node):
                motif_avoidant_count += len(attr)
        if (
            motif_avoidant_count != 0
        ):  # ignore source SCCs with motif avoidant attractors
            exist_maa = True

    # delete the implicit parameters from the node subspaces and the edge motifs
    for node_id in scc_sd.node_ids():
        for implicit in implicit_parameters:
            cast(dict[str, int], scc_sd.dag.nodes[node_id]["space"]).pop(implicit, None)

    for x, y in cast(Iterable[tuple[int, int]], scc_sd.dag.edges):
        for implicit in implicit_parameters:
            cast(dict[str, int], scc_sd.dag.edges[x, y]["motif"]).pop(implicit, None)

    return scc_sd, exist_maa


def attach_scc_sd(
    sd: balm.SuccessionDiagram.SuccessionDiagram,
    scc_sd: balm.SuccessionDiagram.SuccessionDiagram,
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
            int, list(scc_sd.dag.predecessors(scc_node_id))[0]  # type: ignore
        )  # get the first parent
        assert scc_parent_id < scc_node_id

        if scc_parent_id == 0:  # the parent is the root of scc_sd
            parent_id = branch
        else:
            parent_id = size_before_attach + scc_parent_id - 1

        motif = scc_sd.edge_stable_motif(scc_parent_id, scc_node_id)
        motif.update(cast(dict[str, int], sd.dag.nodes[branch]["space"]))

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
            motif.update(cast(dict[str, int], sd.dag.nodes[branch]["space"]))

            child_id = sd._ensure_node(parent_id, motif)  # type: ignore
            assert child_id == size_before_attach + scc_child_id - 1

        # if the node had any child node, consider it expanded.
        if len(scc_child_ids) > 0:
            sd.dag.nodes[parent_id]["expanded"] = True

    return next_branches
