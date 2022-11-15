"""
Modification of Trappist/trappist.py that is easier to invoke directly
as a function, instead of calling it as a separate program.

Also, it allows us to call trappist with restriction to a particular
subspace of the model.
"""

import argparse
import json
import os
import subprocess
import sys
import time
import tempfile
import xml.etree.ElementTree as etree
from typing import Generator, IO, List

from sys import setrecursionlimit

import networkx as nx

from bnet import read_bnet

version = "0.4.1"

setrecursionlimit(2048)


def pnml_to_asp(name: str) -> str:
    """Convert a PNML id to an ASP variable."""
    # TODO handle non-accetable chars
    if name.startswith("-"):
        return "n" + name[1:]
    return "p" + name

def pnml_to_name(name: str) -> str:
    if name.startswith("-"):
        return name[1:]
    return name

def write_asp(petri_net: nx.DiGraph, asp_file: IO, computation: str, time_reversal: str, subspace, source_nodes: List[str]):
    """Write the ASP program for the conflict-free siphons of petri_net."""
    places = []
    free_places = []
    for node, kind in petri_net.nodes(data="kind"):
        if kind == "place":
            places.append(node)
            if not pnml_to_name(node) in subspace:
                free_places.append(node)

            print("{", pnml_to_asp(node), "}.", file=asp_file, sep="")
            if not node.startswith("-"):
                print(
                    f":- {pnml_to_asp(node)}, {pnml_to_asp('-' + node)}.", file=asp_file
                )  # conflict-freeness
        else:  # it's a transition, apply siphon (if one succ is true, one pred must be true)
            if not time_reversal: # compute siphons
                preds = list(petri_net.predecessors(node))
                or_preds = "; ".join(map(pnml_to_asp, preds))
                for succ in petri_net.successors(node):
                    if succ not in preds:  # optimize obvious tautologies
                        print(f"{or_preds} :- {pnml_to_asp(succ)}.", file=asp_file)
            else: # compute traps
                succs = list(petri_net.successors(node))
                or_succs = "; ".join(map(pnml_to_asp, succs))

                for pred in petri_net.predecessors(node):
                    if pred not in succs:  # optimize obvious tautologies
                        print(f"{or_succs} :- {pnml_to_asp(pred)}.", file=asp_file)

    if computation == "max":
        #print(" ".join(free_places))
        max_condition = "; ".join(pnml_to_asp(node) for node in free_places)
        print(
            f"{max_condition}.", file=asp_file
        )

        """More constraints for source nodes"""
        for node in source_nodes:
            if not node in subspace:
                pA = pnml_to_asp(node)
                nA = pnml_to_asp("-" + node)

                print(f"{pA} ; {nA}."
                    , file=asp_file
                )

    for fixed_variable in subspace:
        if subspace[fixed_variable] == "1":
            print(f"n{fixed_variable}.", file=asp_file)            
        if subspace[fixed_variable] == "0":
            print(f"p{fixed_variable}.", file=asp_file)


    if computation == "fix":
        for node in places:
            if not node.startswith("-"):
                print(
                    f"{pnml_to_asp(node)} ; {pnml_to_asp('-' + node)}.", file=asp_file
                )


def solve_asp(asp_filename: str, max_output: int, time_limit: int, computation: str) -> str:
    """Run an ASP solver on program asp_file and get the solutions."""
    dom_mod = "--dom-mod=3, 16" # for min. trap spaces and fixed points

    if computation == "max":
        dom_mod = "--dom-mod=5, 16" # for max. trap spaces
        
    result = subprocess.run(
        [
            "clingo",
            str(max_output),
            "--heuristic=Domain",
            "--enum-mod=domRec",
            dom_mod,
            "--outf=2",  # json output
            f"--time-limit={time_limit}",
            asp_filename,
        ],
        capture_output=True,
        text=True,
    )

    # https://www.mat.unical.it/aspcomp2013/files/aspoutput.txt
    # 30: SAT, all enumerated, optima found, 10 stopped by max, 20 query is false
    if result.returncode != 30 and result.returncode != 10 and result.returncode != 20:
        print(f"Return code from clingo: {result.returncode}")
        result.check_returncode()  # will raise CalledProcessError

    if result.returncode == 20:
        return "UNSATISFIABLE"

    return result.stdout


def solution_to_bool(places: List[str], sol: List[str]) -> List[str]:
    """Convert a list of present places in sol, to a tri-valued vector."""
    return [place_in_sol(sol, p) for p in places]


def place_in_sol(sol: List[str], place: str) -> str:
    """Return 0/1/- if place is absent, present or does not appear in sol.
    Remember that being in the siphon means staying empty, so the opposite value is the one fixed.
    """
    if "p" + place in sol:
        return "0"
    if "n" + place in sol:
        return "1"
    return "-"


def get_solutions(
    asp_output: str, places: List[str]
) -> Generator[List[str], None, None]:
    """Display the ASP output back as trap spaces."""
    solutions = json.loads(asp_output)
    yield from (
        solution_to_bool(places, sol["Value"])
        for sol in solutions["Call"][0]["Witnesses"]
    )


def get_asp_output(
    petri_net: nx.DiGraph, max_output: int, time_limit: int, computation: str, time_reversal, subspace, source_nodes: List[str]
) -> str:
    """Generate and solve ASP file."""
    (fd, tmpname) = tempfile.mkstemp(suffix=".lp", text=True)
    with open(tmpname, "wt") as asp_file:
        write_asp(petri_net, asp_file, computation, time_reversal, subspace, source_nodes)
    solutions = solve_asp(tmpname, max_output, time_limit, computation)

    os.close(fd)
    os.unlink(tmpname)
    return solutions


def compute_trap_spaces(
    petri_net: nx.DiGraph,
    source_nodes: List[str] = [],
    max_output: int = 0,
    time_limit: int = 0,
    computation: str = "min",    
    time_reversal: bool = False,
    subspace = {},
) -> List[str]:
    start = time.process_time()
    
    places = []
    for node, kind in petri_net.nodes(data="kind"):
        if kind == "place" and not node.startswith("-"):
            places.append(node)
    
    computed_object = "min. trap spaces"
    if computation == "max":
        computed_object = "max. trap spaces"
    elif computation == "min":
        computed_object = "min. trap spaces"
    elif computation == "fix":
        computed_object = "fixed points"
    else:
        raise ValueError("Support computing only max. trap spaces, min. trap spaces, and fixed points")

    solutions_output = get_asp_output(petri_net, max_output, time_limit, computation, time_reversal, subspace, source_nodes)

    if solutions_output == "UNSATISFIABLE":
        return []            
    else:
        solutions = get_solutions(solutions_output, places)
        result = []
        for solution in solutions:
            space = {}
            for (place,value) in zip(places, solution):
                if value == "0":
                    space[place] = "0"
                if value == "1":
                    space[place] = "1"
            result.append(space)

        return result
