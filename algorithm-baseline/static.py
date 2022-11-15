import networkx as nx
import matplotlib.pyplot as plt
import random # importing the random module

from pyeda.boolalg import boolfunc
from pyeda.boolalg.bdd import bddvar, expr2bdd
from pyeda.boolalg.expr import expr

from typing import List, Set

"""
A python package for approximating minimum feedback vertex sets
"""
from FVSpython3 import FVS as FVS

def find_minimum_NFVS(network):
    nodes = []
    source_nodes = []
    INx = {}

    for variable in network.variables():
        var_name = network.get_variable_name(variable)
        function = network.get_update_function(variable)

        nodes.append(var_name)
        if function.strip() == var_name:
            source_nodes.append(var_name)

        fx = expr(function.replace("!", "~"))

        INx[var_name] = fx.support # list of nodes appearing in Boolean function fx
            
    """Build the unsigned interaction graph"""
    ig = nx.DiGraph()

    for x in nodes:
        ig.add_node(x)

        for y in INx[x]:
            ig.add_edge(str(y), str(x)) 

        """
        TODO: node y may not actually affect node x
        We may need to use BDDs here.
        """

    """Find feedback vertex set"""
    U = FVS.FVS(ig)

    U_neg = U # TODO: compute a minimum negative feedback vertex set

    return U_neg


