from typing import IO

import networkx as nx
import re

from pyeda.boolalg import boolfunc
from pyeda.boolalg.bdd import bddvar, expr2bdd
from pyeda.boolalg.expr import expr

def space_to_string(network, space):
    result = ""
    for var in network.variables():
        name = network.get_variable_name(var)
        if name in space:
            result += space[name]
        else:
            result += "-"
    return result

def space_to_aeon_set(stg, space):
    """
        Converts a space (dictionary of variable values) to a symbolic BDD represented set.

        TODO: AEON Rust has a faster function for this, but it is not exported to Python yet.
        Once that is done, this should be much easier to implement.
    """
    result = stg.unit_colored_vertices()
    for var_name in space:
        value = space[var_name]
        if value == "1":
            result = result.intersect(stg.fix_variable(var_name, True))
        if value == "0":
            result = result.intersect(stg.fix_variable(var_name, False))
    return result


def add_edges(net: nx.DiGraph, tname: str, things: boolfunc, source: str):
    """Add all edges from things to tname and back in net except for source."""
    if source.startswith("-"):
        nsource = source[1:]
    else:
        nsource = "-" + source

    for i, t in enumerate(things.satisfy_all()):
        name = f"{tname}_{i}"        
        net.add_node(name, kind="transition")
        for p, v in t.items():
            if v == 0:
                pname = "-" + str(p)
            else:
                pname = str(p)
            net.add_edge(pname, name)
            if pname == source:
                net.add_edge(name, nsource)
            else:
                net.add_edge(name, pname)

def aeon_to_petri_net(network):
    """
        Takes an AEON representation of a Boolean network and converts it to DiGraph
        for processing by Trappist.

        TODO: At the moment, this uses pyeda BDDs for conversion, but AEON BDDs should be
        simpler once necessary APIs are in place.
    """
    net = nx.DiGraph()

    for variable in network.variables():
        var_name = network.get_variable_name(variable)
        function = network.get_update_function(variable)

        net.add_node(var_name, kind="place")
        net.add_node("-" + var_name, kind="place")
        vx = bddvar(var_name)
        fx = expr(function.replace("!", "~"))
        fx = expr2bdd(fx)
        activate = fx & ~vx
        inactivate = ~fx & vx

        add_edges(net, f"tp_{var_name}", activate, "-" + var_name)
        add_edges(net, f"tn_{var_name}", inactivate, var_name)

    return net
