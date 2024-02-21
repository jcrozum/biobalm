import networkx as nx
import numpy as np
import scipy as sp


def power_law_graph_generator(n_nodes: int, power: float = 3.0, seed: int = 0):
    """
    Function to generate a random graph with power-law out-degree distribution
    and a Poisson in-degree distribution.
    """
    rng = np.random.default_rng(seed)

    A = 1 / sp.special.zeta(power)  # normalization factor

    ks = []
    for i in range(n_nodes):
        rand = 1 - rng.random()  # (0.0, 1.0]
        k = 0
        val = 0
        while val < rand:
            k += 1
            val += A / k**power
            # reset and restart if the degree gets bigger than the number of nodes
            # this cuts off the tail of the power law distribution
            if k > n_nodes:
                rand = 1 - rng.random()  # (0.0, 1.0]
                k = 0
                val = 0
        ks.append(k)

    # print(ks)
    # print(np.average(ks))

    G = nx.DiGraph()
    while np.sum(ks) > 0:  # type: ignore
        source, sink = rng.integers(0, n_nodes, 2)
        # if all the out degree is filled, find a new source
        while ks[source] == 0:
            source = rng.integers(0, n_nodes)
        # if the edge already exists, find a new sink
        while G.has_edge("n" + str(source), "n" + str(sink)):
            sink = rng.integers(0, n_nodes)

        G.add_edge("n" + str(source), "n" + str(sink))
        ks[source] = ks[source] - 1

    # nx.draw(G)
    # plt.show()

    # degrees = list(d for n, d in G.degree())
    # print(np.average(degrees))

    # G2 = G.copy()
    # removedNodes = []
    # largestSCC = max(nx.strongly_connected_components(G2), key = len)
    # for node in G2.nodes():
    #     if(node not in largestSCC):
    #         removedNodes.append(node)
    # G2.remove_nodes_from(removedNodes)

    return G


def add_negative_edges(G: nx.DiGraph, pNeg: float = 0.25, seed: int = 0):
    """
    Function to add negative edges to graph G
    inputs:
        G = Graph to add negative edges to
        pNeg = Probability that an edge is negative
        seed = Random number generator seed
    """
    rng = np.random.default_rng(seed)
    for E in G.edges():
        isNeg = rng.random()
        if isNeg < pNeg:
            G.edges()[E]["negative"] = True
        else:
            G.edges()[E]["negative"] = False


def generate_ncf_rule(
    G: nx.DiGraph, node: int, bias: float = 0.5, seed: int = 0
) -> str:
    """
    Function to write a nested canalyzing function for a node as a read-once function (RoF)
    inputs:
        G = Graph
        node = Node to generate rule for
        rng = Random number generator
        bias = Probability that the function input is sufficient (or sufficient inhibitor)
               of the output at each canalyzing depth. Also determines the ratio of 1 for
               constant nodes.
    output:
        String of nested canalyzing function for given node
    """
    rng = np.random.default_rng(seed)

    inputs = [e[0] for e in G.in_edges(node)]
    negative = [G.get_edge_data(x, node)["negative"] for x in inputs]
    n = len(inputs)
    if n == 0:
        rand = rng.random()
        val = "false"
        if rand < bias:
            val = "true"
        return str(node) + ", " + val
    order = list(rng.choice(inputs, n, replace=False))
    funct = str(node) + ", "
    for i in range(len(order)):
        rand = rng.random()
        if rand < bias:
            if negative[inputs.index(order[i])]:
                funct += "!" + str(order[i]) + " & ("
            else:
                funct += str(order[i]) + " | ("
        else:
            if negative[inputs.index(order[i])]:
                funct += "!" + str(order[i]) + " | ("
            else:
                funct += str(order[i]) + " & ("
    funct = funct.rstrip("& (|")
    funct += ")" * (len(inputs) - 1)
    return funct
