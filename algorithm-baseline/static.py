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

def get_source_nodes(network):
    source_nodes = []

    for variable in network.variables():
        var_name = network.get_variable_name(variable)
        function = network.get_update_function(variable)

        if function.strip() == var_name:
            source_nodes.append(var_name)

    return source_nodes

def find_minimum_NFVS(network):
    nodes = []
    INx = {}
    funs_bdd = {}
    vars_bdd = {}

    for variable in network.variables():
        var_name = network.get_variable_name(variable)
        function = network.get_update_function(variable)

        nodes.append(var_name)

        fx = expr(function.replace("!", "~"))

        INx[var_name] = fx.support # list of nodes appearing in Boolean function fx

        vx = bddvar(var_name)
        fx = expr2bdd(expr(fx))

        vars_bdd[var_name] = vx
        funs_bdd[var_name] = fx
            
    """Build the unsigned and signed interaction graphs"""
    u_ig = nx.DiGraph()
    s_ig = nx.DiGraph()

    for x in nodes:
        u_ig.add_node(x)
        s_ig.add_node(x)

        fx = funs_bdd[x]

        for y in INx[x]:
            is_actual_arc = False

            vy = vars_bdd[str(y)]
            
            fx_res_vy_0 = fx.restrict({vy: 0})
            fx_res_vy_1 = fx.restrict({vy: 1})

            if (~fx_res_vy_0 & fx_res_vy_1).satisfy_one():
                # a positive arc with weight = 1
                s_ig.add_edge(str(y), str(x), weight=1)
                is_actual_arc = True
                print(str(y) + "<1>" + str(x))

            if (fx_res_vy_0 & ~fx_res_vy_1).satisfy_one():
                # a negative arc with weight = -1
                s_ig.add_edge(str(y), str(x), weight=-1)
                is_actual_arc = True
                print(str(y) + "<-1>" + str(x))

            if is_actual_arc == True:
                u_ig.add_edge(str(y), str(x))
                #print(str(y) + "<>" + str(x))


    """Find feedback vertex set"""
    U = FVS.FVS(u_ig)

    print(U)

    U_neg = U # TODO: compute a minimum negative feedback vertex set

    return U_neg


