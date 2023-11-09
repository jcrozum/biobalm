import sys
sys.path.append(".")

import networkx as nx
import itertools as it

from biodivine_aeon import BooleanNetwork
from nfvsmotifs.interaction_graph_utils import infer_signed_interaction_graph
from nfvsmotifs.space_utils import percolate_space, percolate_network
from nfvsmotifs.SuccessionDiagram import SuccessionDiagram
from networkx import DiGraph
from nfvsmotifs.petri_net_translation import extract_variable_names, network_to_petrinet

DEBUG = False

def expand_source_SCCs(sd: SuccessionDiagram) -> bool:
    """
    1. percolate
    2. find source nodes, fix combinations and then percolate
    3. at every level, find all source SCCs, expand them
    4. when there are no more SCCs, expand in a usual way.

    TODO: it will be rare especially for empirical models,
    but there can be source SCCs even before fixing source node combinations.
    It could be useful to find them first, calculate scc_sd and store them somewhere,
    rather than calculating the same thing for every source node combination
    """

    root = sd.root()

    current_level = [root]
    next_level: list[int] = []
    final_level: list[int] = [] # from here there are no more source SCCs

    # percolate constant nodes
    perc_space, _ = percolate_space(sd.network, {}, strict_percolation=False)
    sd.G.nodes[root]["space"] = perc_space

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

            next_level.append(sd._ensure_node(root, sub_space))
        
        sd.G.nodes[root]["expanded"] = True
        sd.G.nodes[root]["attractors"] = [] # no need to look for attractors here
        current_level = next_level
        next_level = []

    # each level consists of one round of fixing all source SCCs
    while len(current_level) > 0:
        if DEBUG:
            print(f"{current_level=}")
        for node_id in current_level:
            sub_space = sd.G.nodes[node_id]["space"]
            
            # find source SCCs
            clean_bnet, clean_bn = perc_and_remove_constants_from_bn(perc_bn, sub_space)
            source_scc_list = find_source_SCCs(clean_bn)

            if DEBUG:
                print(f"{source_scc_list=}")

            if len(source_scc_list) == 0: # no more source SCCs
                final_level.append(node_id)
                continue
            
            sd.G.nodes[node_id]["attractors"] = [] # no need to look for attractors here
            current_branches = [node_id] # this is where the scc_sd should be "attached"
            next_branches = []

            while len(source_scc_list) > 0:
                source_scc = source_scc_list.pop(0)
                scc_sd = find_scc_sd(clean_bnet, source_scc)

                # check for motif avoidant attractors
                motif_avoidant_count = 0
                for node in scc_sd.node_ids():
                    attr = scc_sd.node_attractor_seeds(node, compute=True)
                    if not scc_sd.node_is_minimal(node):
                        motif_avoidant_count += len(attr)
                if motif_avoidant_count != 0: # ignore source SCCs with motif avoidant attractors
                    #TODO: somehow skip this calculation when this source SCC appears again later.
                    continue
                
                # add scc_sd to every single branch in the original sd
                for branch in current_branches:
                    new_branches = attach_scc_sd(sd, scc_sd, branch)
                    next_branches.extend(new_branches)

                current_branches = next_branches
                next_branches = []
            
            if current_branches == [node_id]:
                # This happens when there were only source SCCs with motif avoidant attractors
                final_level.append(node_id)
            else:
                # once all the source SCCs are expanded, the final branches become the next level
                next_level.extend(current_branches)
        
        current_level = next_level
        next_level = []

    # Finally, expand the sd in a usual way.
    # TODO: motif avoidance that is already calculated in the source SCCs should be used here somehow.
    if DEBUG:
        print(f"{final_level=}")
    
    for node_id in final_level:
        sd.G.nodes[node_id]["attractors"] = None # check attractors from here
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


def perc_and_remove_constants_from_bn(bn:BooleanNetwork, space: dict[str, int]) -> tuple[str, BooleanNetwork]:
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

    perc_space, _ = percolate_space(bn, space, strict_percolation=False)
    perc_bn = percolate_network(bn, perc_space)

    perc_bnet = perc_bn.to_bnet()

    # remove constant nodes from the bnet
    clean_bnet = ""
    for line in perc_bnet.split("\n"):
        if line.strip() == "":
            continue
        node = line.split(",")[0]
        function = line.split(",")[1].strip()
        if function == "("+node+" | !"+node+")":
            continue
        elif function == "("+node+" & !"+node+")":
            continue
        else:
            clean_bnet += line + "\n"

    clean_bn = BooleanNetwork.from_bnet(clean_bnet)

    return clean_bnet, clean_bn


def find_source_SCCs(bn:BooleanNetwork) -> list[list[str]]:
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
    scc_list = [sorted(x) for x in nx.strongly_connected_components(di)]

    # find source SCCs
    source_scc_list = []
    for scc in scc_list:
        source = True
        for node in scc:
            for reg in di.predecessors(node):
                if reg not in scc:
                    source = False
                    break
            if not source:
                break
        if source:
            source_scc_list.append(scc)

    return sorted(source_scc_list)


def find_scc_sd(bnet:str, source_scc:list[str]) -> SuccessionDiagram:
    """
    TODO: better way that does not use bnet but rather bn directly or petri_net directly
    to find the scc_sd would be useful.
    TODO: somehow remove implicit parameters elegantly
    """

    # get rules for only the source SCC
    scc_bnet = "targets,factors\n"
    for line in bnet.split("\n"):
        for node in source_scc:
            if line.startswith(node+","):
                scc_bnet += line + "\n"

    scc_bn = BooleanNetwork.from_bnet(scc_bnet)

    if DEBUG:
        print("scc_bnet\n",scc_bn.to_bnet())

    # scc_bn = scc_bn.infer_regulatory_graph()

    # quick way to ignore implicit parameters coming from fake edges into source SCCs
    # e.g. in a rule such as A* = A | A & B, B may appear in the scc_bnet implicitly but is meaningless
    implicit_parameters = []
    for var in scc_bn.implicit_parameters():
        scc_bn.set_update_function(var, "true")
        implicit_parameters.append(scc_bn.get_variable_name(var))
    scc_bn = scc_bn.infer_regulatory_graph()

    # Compute the succession diagram.
    scc_sd = SuccessionDiagram(scc_bn)
    fully_expanded = scc_sd.expand_bfs()
    assert fully_expanded

    # delete the implicit parameters from the node subspaces
    # unfortunately the edges still carry the implicit parameters, but they are not used for the full sd anyways
    for node_id in scc_sd.node_ids():
        for implicit in implicit_parameters:
            scc_sd.G.nodes[node_id]["space"].pop(implicit, None)

    return scc_sd

def attach_scc_sd(sd:SuccessionDiagram, scc_sd:SuccessionDiagram, branch:int) -> list[int]:
    """
    Attach scc_sd to the given branch point of the sd.
    Returns the new branching points.
    Make sure to only use scc_sd that does not have motif avoidant attractors,
    as the nodes will be assumed to not have any attractors
    """
    if len(scc_sd) == 1:
        return [branch]

    next_branches = []
    current_size = sd.G.number_of_nodes()
    # first add all the nodes using their first parent
    for scc_node_id in scc_sd.node_ids():
        if scc_node_id == 0: # no need to add the root of scc_sd
            continue

        scc_parent_id = list(scc_sd.G.predecessors(scc_node_id))[0] # get the first parent
        assert scc_parent_id < scc_node_id

        if scc_parent_id == 0: # the parent is the root of scc_sd
            parent_id = branch
        else:
            parent_id = current_size + scc_parent_id - 1

        scc_sub_space = {}
        scc_sub_space.update(sd.G.nodes[branch]["space"])
        scc_sub_space.update(scc_sd.G.nodes[scc_node_id]["space"])

        child_id = sd._ensure_node(parent_id, scc_sub_space)
        sd.G.nodes[child_id]["attractors"] = [] # no need to check for attractors in these nodes

        if scc_sd.node_is_minimal(scc_node_id) and child_id not in next_branches:
            next_branches.append(child_id)

    # now add all the missing edges
    for scc_node_id in scc_sd.node_ids():
        if scc_node_id == 0:
            node_id = branch
        else:
            node_id = current_size + scc_node_id - 1
        for scc_child_id in scc_sd.G.successors(scc_node_id):
            scc_sub_space = {}
            scc_sub_space.update(sd.G.nodes[branch]["space"])
            scc_sub_space.update(scc_sd.G.nodes[scc_child_id]["space"])
            
            child_id = sd._ensure_node(node_id, scc_sub_space)
            assert child_id == current_size + scc_child_id - 1

        sd.G.nodes[node_id]["expanded"] = True

    return next_branches