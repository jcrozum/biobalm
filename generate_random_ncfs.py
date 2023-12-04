import os

from rbn_generators import *

N_NODES_LIST = [10, 20, 40, 80, 160, 320]
GRAPHS_PER_SIZE = 100
POWER = 2.5  # power for the out-degree, diverges at 2.0
P_NEG = 0.25  # ratio of negative edges
DIRECTORY = "models/random_ncf2"

for N_NODES in N_NODES_LIST:
    for i in range(GRAPHS_PER_SIZE):
        G = power_law_graph_generator(N_NODES, POWER, seed=i)

        add_negative_edges(G, P_NEG, seed=i)

        rules = ""
        for node in G.nodes():
            rules += generate_ncf_rule(G, node, bias=0.5, seed=i) + "\n"
        writePath = DIRECTORY + "/n" + str(N_NODES) + "_" + str(i) + ".bnet"
        isExist = os.path.exists(DIRECTORY)
        if not isExist:
            os.makedirs(DIRECTORY)

        with open(writePath, "w") as f:
            f.write("targets, factors\n")
            f.write(rules)
