import sys

sys.path.append("..")

from biodivine_aeon import BooleanNetwork

import nfvsmotifs
import networkx as nx
from nfvsmotifs.interaction_graph_utils import infer_signed_interaction_graph, _digraph_to_regulatory_graph
from nfvsmotifs.space_utils import percolate_space, percolate_network
from nfvsmotifs.SuccessionDiagram import SuccessionDiagram

# This is unfortunately necessary for PyEDA Boolean expression parser (for now).
sys.setrecursionlimit(150000)

NODE_LIMIT = 1_000
DEPTH_LIMIT = 100

LOG_LOCATION = "SCC_analysis_" + sys.argv[1].split("/")[-2] + ".csv"

bn = BooleanNetwork.from_file(sys.argv[1])

perc_space, conflict = percolate_space(bn, {}, strict_percolation=False)
perc_bn = percolate_network(bn, perc_space)

bnet = perc_bn.to_bnet()

perc_bnet = ""
for line in bnet.split("\n"):
    if line.strip() == "":
        continue
    node = line.split(",")[0]
    function = line.split(",")[1].strip()
    if function == "("+node+" | !"+node+")":
        continue
    elif function == "("+node+" & !"+node+")":
        continue
    else:
        perc_bnet += line + "\n"

perc_bn = BooleanNetwork.from_bnet(perc_bnet)

di = infer_signed_interaction_graph(perc_bn)

scc_list = [x for x in nx.strongly_connected_components(di)]

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

print(source_scc_list)

print(perc_bnet)

est_size = 1
est_depth = 0
est_min = 1

for source_scc in source_scc_list:
    print("source_scc:", source_scc)
    scc_bnet = "targets,factors\n"
    for line in perc_bnet.split("\n"):
        for node in source_scc:
            if line.startswith(node+","):
                scc_bnet += line + "\n"

    scc_bn = BooleanNetwork.from_bnet(scc_bnet)

    print(scc_bn.to_bnet())
    scc_bn = scc_bn.infer_regulatory_graph()

    for var in scc_bn.implicit_parameters():
        scc_bn.set_update_function(var, "true")

    scc_bn = scc_bn.infer_regulatory_graph()


    # Compute the succession diagram.
    sd = SuccessionDiagram(scc_bn)
    fully_expanded = sd.expand_bfs(bfs_level_limit=DEPTH_LIMIT, size_limit=NODE_LIMIT)
    assert fully_expanded

    print(f"SD size, {len(sd)}")
    print(f"SD depth, {sd.depth()} ")
    print(f"minimal trapspaces, {len(sd.minimal_trap_spaces())}")

    est_size *= len(sd)
    est_depth += sd.depth()
    est_min *= len(sd.minimal_trap_spaces())

print(f"{est_size=}")
print(f"{est_depth=}")
print(f"{est_min=}")

log = open(LOG_LOCATION, "a")
log.write(sys.argv[1].split("/")[-1] + ",")  # model name
log.write(str(bn.num_vars()) + ",")  # network size
log.write(str(len(source_scc_list)) + ",")  # number of source_scc
log.write(str(len(scc_list)) + ",") # number of scc
log.write(str(est_size) + ",") # estimated SD size
log.write(str(est_depth) + ",") # estimated SD depth
log.write(str(est_min) + "\n") # estimated number of attractors
