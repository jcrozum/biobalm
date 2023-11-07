import sys

sys.path.append(".")

from biodivine_aeon import BooleanNetwork

import nfvsmotifs
import networkx as nx
import itertools as it
from nfvsmotifs.interaction_graph_utils import infer_signed_interaction_graph, _digraph_to_regulatory_graph
from nfvsmotifs.space_utils import percolate_space, percolate_network
from nfvsmotifs.SuccessionDiagram import SuccessionDiagram
from networkx import DiGraph
from nfvsmotifs.petri_net_translation import extract_variable_names, network_to_petrinet

def expand_source_SCCs(sd: SuccessionDiagram, size_limit: int | None = None) -> bool:
    """
    See `SuccessionDiagram.expand_minimal_spaces` for documentation.
    """

    root = sd.root()

    seen = set([root])

    # percolate constant nodes
    perc_space, _ = percolate_space(sd.network, {}, strict_percolation=False)
    perc_bn = percolate_network(sd.network, perc_space)

    # find source nodes
    source_nodes = find_source_nodes(perc_bn)

    # get the initial subspaces that constitute the first depth
    initial_sub_spaces:list[dict[str,int]] = []
    if len(source_nodes) == 0:
        initial_sub_spaces.append(perc_space)
    else:
        # get source nodes combinations
        bin_values_list = it.product(range(2), repeat=len(source_nodes)) 
        for bin_values in bin_values_list:
            source_comb = dict(zip(source_nodes, bin_values))

            sub_space = source_comb
            sub_space.update(perc_space)
            initial_sub_spaces.append(sub_space)


def perc_and_remove_constants_from_bn(bn:BooleanNetwork, space: dict[str, int]) -> BooleanNetwork:
    """
    Take a BooleanNetwork and percolate given space.
    Then remove constant nodes with rules
    A* = (A | !A)
    B* = (B & !B)
    and return a clean BooleanNetwork with no constant nodes.

    TODO: if possible, delete constant nodes from the BooleanNetwork directly,
    without converting into bnet and back.
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

    return clean_bn


def find_source_nodes(network: BooleanNetwork | DiGraph) -> list[str]:
    
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

# di = infer_signed_interaction_graph(perc_bn)

# scc_list = [x for x in nx.strongly_connected_components(di)]

# source_scc_list = []
# for scc in scc_list:
#     source = True
#     for node in scc:
#         for reg in di.predecessors(node):
#             if reg not in scc:
#                 source = False
#                 break
#         if not source:
#             break
#     if source:
#         source_scc_list.append(scc)

# print(source_scc_list)

# print(perc_bnet)

# est_size = 1
# est_depth = 0
# est_min = 1

# for source_scc in source_scc_list:
#     print("source_scc:", source_scc)
#     scc_bnet = "targets,factors\n"
#     for line in perc_bnet.split("\n"):
#         for node in source_scc:
#             if line.startswith(node+","):
#                 scc_bnet += line + "\n"

#     scc_bn = BooleanNetwork.from_bnet(scc_bnet)

#     print(scc_bn.to_bnet())
#     scc_bn = scc_bn.infer_regulatory_graph()

#     for var in scc_bn.implicit_parameters():
#         scc_bn.set_update_function(var, "true")

#     scc_bn = scc_bn.infer_regulatory_graph()


#     # Compute the succession diagram.
#     sd = SuccessionDiagram(scc_bn)
#     fully_expanded = sd.expand_bfs(bfs_level_limit=DEPTH_LIMIT, size_limit=NODE_LIMIT)
#     assert fully_expanded

#     print(f"SD size, {len(sd)}")
#     print(f"SD depth, {sd.depth()} ")
#     print(f"minimal trapspaces, {len(sd.minimal_trap_spaces())}")

#     est_size *= len(sd)
#     est_depth += sd.depth()
#     est_min *= len(sd.minimal_trap_spaces())

# print(f"{est_size=}")
# print(f"{est_depth=}")
# print(f"{est_min=}")

# log = open(LOG_LOCATION, "a")
# log.write(sys.argv[1].split("/")[-1] + ",")  # model name
# log.write(str(bn.num_vars()) + ",")  # network size
# log.write(str(len(source_scc_list)) + ",")  # number of source_scc
# log.write(str(len(scc_list)) + ",") # number of scc
# log.write(str(est_size) + ",") # estimated SD size
# log.write(str(est_depth) + ",") # estimated SD depth
# log.write(str(est_min) + "\n") # estimated number of attractors
