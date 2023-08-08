from rbn_generators import *
import os

N_GRAPHS = 20
N_NODES = 300
POWER = 2.5     # power for the out-degree, diverges at 2.0
P_NEG = 0.25    # ratio of negative edges
DIRECTORY = 'random_ncf'

for i in range(N_GRAPHS):
    G = power_law_graph_generator(N_NODES, POWER, seed = i)

    add_negative_edges(G, P_NEG, seed = i)

    rules = ""
    for node in G.nodes():
        rules += generate_ncf_rule(G, node, bias = 0.5, seed = i)+'\n'
    writePath = DIRECTORY+'/n'+ str(N_NODES) +'_'+str(i)+'.bnet'
    isExist = os.path.exists(DIRECTORY)
    if not isExist:
        os.makedirs(DIRECTORY)       
    f = open(writePath,'w')
    f.write('targets, factors\n')
    f.write(rules)
    f.close()