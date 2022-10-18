"""Compute trap spaces of a Petri-net encoded Boolean model.

Copyright (C) 2022 Sylvain.Soliman@inria.fr and giang.trinh91@gmail.com

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
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

def read_pnml(fileobj: IO) -> nx.DiGraph:
    """Parse the given file."""
    root = etree.parse(fileobj).getroot()
    if root.tag != "pnml":
        raise ValueError("Currently limited to parsing PNML files")
    net = nx.DiGraph()

    for place in root.findall("./net/place"):
        net.add_node(
            place.get("id"), kind="place"
        )

    for transition in root.findall("./net/transition"):
        net.add_node(transition.get("id"), kind="transition")

    for arc in root.findall("./net/arc"):
        net.add_edge(arc.get("source"), arc.get("target"))

    return net


def pnml_to_asp(name: str) -> str:
    """Convert a PNML id to an ASP variable."""
    # TODO handle non-accetable chars
    if name.startswith("-"):
        return "n" + name[1:]
    return "p" + name


def write_asp(petri_net: nx.DiGraph, asp_file: IO, computation: str):
    """Write the ASP program for the conflict-free siphons of petri_net."""
    places = []
    for node, kind in petri_net.nodes(data="kind"):
        if kind == "place":
            places.append(node)
            print("{", pnml_to_asp(node), "}.", file=asp_file, sep="")
            if not node.startswith("-"):
                print(
                    f":- {pnml_to_asp(node)}, {pnml_to_asp('-' + node)}.", file=asp_file
                )  # conflict-freeness
        else:  # it's a transition, apply siphon (if one succ is true, one pred must be true)
            preds = list(petri_net.predecessors(node))
            or_preds = "; ".join(map(pnml_to_asp, preds))
            for succ in petri_net.successors(node):
                if succ not in preds:  # optimize obvious tautologies
                    print(f"{or_preds} :- {pnml_to_asp(succ)}.", file=asp_file)
                    
    if computation == "max":
        max_condition = "; ".join(pnml_to_asp(node) for node in places)
        print(
            f"{max_condition}.", file=asp_file
        )

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
    petri_net: nx.DiGraph, max_output: int, time_limit: int, computation: str
) -> str:
    """Generate and solve ASP file."""
    (_, tmpname) = tempfile.mkstemp(suffix=".lp", text=True)
    with open(tmpname, "wt") as asp_file:
        write_asp(petri_net, asp_file, computation)
    solutions = solve_asp(tmpname, max_output, time_limit, computation)
    os.unlink(tmpname)
    return solutions


def compute_trap_spaces(
    infile: IO,
    display: bool = False,
    max_output: int = 0,
    time_limit: int = 0,
    computation: str = "min",
    method: str = "asp",
) -> Generator[List[str], None, None]:
    """Do the minimal trap space computation on input file infile."""
    toclose = False
    if isinstance(infile, str):
        infile = open(infile, "r", encoding="utf-8")
        toclose = True

    start = time.process_time()
    if infile.name.endswith(".bnet"):
        petri_net = read_bnet(infile)
    else:
        infile.close()
        raise ValueError("Currently limited to parsing .bnet files")

    if toclose:
        infile.close()

    #print(f"Petri net time {time.process_time() - start:.2f}\n===")
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

    print(f"Compute {computed_object}")

    #print("# nodes = ", len(places))
    if display:
        print(" ".join(places))

    if method == "asp":
        solutions_output = get_asp_output(petri_net, max_output, time_limit, computation)

        if solutions_output == "UNSATISFIABLE":
            print(f"No {computed_object}")
            return
        else:
            solutions = get_solutions(solutions_output, places)

    if display:
        print("\n".join(" ".join(sol) for sol in solutions))
        # print("Total time:", solutions["Time"]["Total"], "s")
        return
    else:
        yield from solutions


def main():
    """Read the Petri-net send the output to ASP and print solution."""
    parser = argparse.ArgumentParser(
        description=" ".join(__doc__.splitlines()[:3]) + " GPLv3"
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s v{version}".format(version=version),
    )
    parser.add_argument(
        "-m",
        "--max",
        type=int,
        default=0,
        help="Maximum number of solutions (0 for all).",
    )
    parser.add_argument(
        "-t",
        "--time",
        type=int,
        default=0,
        help="Maximum number of seconds for search (0 for no-limit).",
    )
    parser.add_argument(
        "-c",
        "--computation",
        choices=["min", "max", "fix"],
        default="min",
        type=str,
        help="Computation option.",
    )
    parser.add_argument(
        "-s",
        "--solver",
        choices=["asp"],
        default="asp",
        type=str,
        help="Solver to compute the conflict-free siphons.",
    )
    parser.add_argument(
        "infile",
        type=argparse.FileType("r", encoding="utf-8"),
        nargs="?",
        default=sys.stdin,
        help="Petri-net (PNML) file",
    )
    args = parser.parse_args()

    try:
        next(compute_trap_spaces(
            args.infile,
            display=True,
            max_output=args.max,
            time_limit=args.time,
            computation=args.computation,
            method=args.solver,
        ))
    except StopIteration:
        pass

if __name__ == "__main__":
    main()

