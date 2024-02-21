import os

import pystablemotifs.random_boolean_networks as rbn

N_NODES_LIST = [10, 20, 40, 80, 160, 320]
GRAPHS_PER_SIZE = 100
K = 2  # in-degree
p = 0.5  # ratio of 1 in the truth table
DIRECTORY = "models/random_nk"

for N_NODES in N_NODES_LIST:
    rules_list = rbn.random_boolean_network_ensemble_kauffman(
        N_NODES, K, p, GRAPHS_PER_SIZE, seed=1000
    )

    for i in range(GRAPHS_PER_SIZE):
        rules = rules_list[i]

        rules = rules.replace(" *=", ",")
        rules = rules.replace(" and ", " & ")
        rules = rules.replace(" or ", " | ")
        rules = rules.replace("not ", "!")
        rules = rules.replace(", 1", ", true")
        rules = rules.replace(", 0", ", false")

        writePath = DIRECTORY + "/n" + str(N_NODES) + "_" + str(i) + ".bnet"
        isExist = os.path.exists(DIRECTORY)
        if not isExist:
            os.makedirs(DIRECTORY)

        with open(writePath, "w") as f:
            f.write("targets, factors\n")
            f.write(rules)
