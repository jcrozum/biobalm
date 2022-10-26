"""
A place for prototyping attractor detection algorithms.

Basic terminology:

 - A space is a dictionary where variable names are its keys and values are strings "0"/"1".
 If a variable is free, it is not present in the dictionary.

 - Percolation refers to the syntactic process of propagating variables fixed in a space using the update
 functions of the specific BN.

"""
from biodivine_aeon import *
from percolation import percolate
from pathlib import Path
from pyeda.boolalg.expr import expr
from trappist import compute_trap_spaces
from conversions import aeon_to_petri_net, space_to_aeon_set
from bnet import read_bnet
from aeon_utils import remove_static_constraints

import sys

def symbolic_attractor_search(stg, candidates):
    while not candidates.is_empty():
        pivot = candidates.pick_vertex()
        bwd = reach_bwd(stg, pivot, candidates)
        fwd = reach_fwd(stg, pivot, candidates)

        scc = fwd.intersect(bwd)
        if fwd.minus(bwd).is_empty():
            print("Attractor", scc)

        candidates = candidates.minus(bwd)


def attractors(stg, petri_net, space):
    space = percolate(stg.network(), space)
    print("Start search in", space)

    universe = space_to_aeon_set(stg, space)
    basins = stg.empty_colored_vertices()

    if len(space) != stg.network().num_vars():        
        max_traps = compute_trap_spaces(petri_net, computation="max", subspace=space)        
        print(f"Found {len(max_traps)} inner trap spaces.")
        for trap in max_traps:
            print(trap)
        for trap in max_traps:
            if trap != space:
                attractors(stg, petri_net, trap)
                basins = basins.union(space_to_aeon_set(stg, trap))

    covered = reach_bwd(stg, basins, universe)

    candidates = universe.minus(covered)
    print("Non-trap candidates:", candidates)
    symbolic_attractor_search(stg, candidates)



if __name__ == "__main__":

    # Load model from the first argument.
    model_path = sys.argv[1]
    if model_path.endswith(".aeon"):
        network = BooleanNetwork.from_aeon(Path(model_path).read_text())
    elif model_path.endswith(".bnet"):
        network = BooleanNetwork.from_bnet(Path(model_path).read_text())
    elif model_path.endswith(".sbml"):
        network = BooleanNetwork.from_sbml(Path(model_path).read_text())
    else:
        print("Unknown model format:", model_path)
        sys.exit()

    # Disable static analysis in AEON (some randomly generated or 
    # pre-processed files might fail this check).
    network = remove_static_constraints(network)
    print(network)

    petri_net = aeon_to_petri_net(network)    
    print(petri_net)

    stg = SymbolicAsyncGraph(network)
    attractors(stg, petri_net, {})
