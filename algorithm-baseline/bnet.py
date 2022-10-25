"""
Modification of Trappist/bnet.py that can accept strings as inputs instead of files.
"""

from typing import IO

import networkx as nx
import re

from pyeda.boolalg import boolfunc
from pyeda.boolalg.bdd import bddvar, expr2bdd
from pyeda.boolalg.expr import expr

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


def read_bnet(string) -> nx.DiGraph:
    """Parse a BoolNet .bnet file and build the corresponding Petri net."""
    net = nx.DiGraph()

    for line in string.split('\n'):        
        if line.startswith("#") or re.match(r'[tT]argets,\s*[Ff]actors', line):
            continue
        try:
            x, fx = line.replace(" ", "").replace("!", "~").split(",", maxsplit=1)
        except ValueError:
            continue
        
        x = x.strip()
        
        net.add_node(x, kind="place")
        net.add_node(
            "-" + x, kind="place"
        )  # convention in PNML files obtained from bnet
        vx = bddvar(x)
        fx = expr2bdd(expr(fx))
        activate = fx & ~vx
        inactivate = ~fx & vx

        add_edges(net, f"tp_{x}", activate, "-" + x)
        add_edges(net, f"tn_{x}", inactivate, x)

    return net
