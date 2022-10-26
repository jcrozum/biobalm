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
from conversions import aeon_to_petri_net, space_to_aeon_set, space_to_string
from bnet import read_bnet
from aeon_utils import remove_static_constraints, has_parameters

import sys

def symbolic_attractor_search(stg, candidates):
    while not candidates.is_empty():
        pivot = candidates.pick_vertex()
        bwd = reach_bwd(stg, pivot, candidates)
        fwd = reach_fwd(stg, pivot, candidates)

        scc = fwd.intersect(bwd)
        if fwd.minus(bwd).is_empty():
            print("Attractor", scc)
            if scc.cardinality() < 100:
                for state in scc.vertices().list_vertices():
                    print(state)

        candidates = candidates.minus(bwd)

def is_subspace(x, y):
    for key in y:
        if (not key in x) or x[key] != y[key]:
            return False
    return True

def attractors(network):
    petri_net = aeon_to_petri_net(network)
    stg = SymbolicAsyncGraph(network)

    def attractors_recursive(space, candidates):
        space_percolated = percolate(network, space)
        if space_percolated != space:
            space = space_percolated
            candidates = candidates.intersect(space_to_aeon_set(stg, space))
            
        universe = space_to_aeon_set(stg, space)
        all_traps = stg.empty_colored_vertices()

        if len(space) != stg.network().num_vars():        
            # Only check for subspaces when the current space still has some free variables remaining.
            # Otherwise the result will always be UNSAT.
            max_traps = compute_trap_spaces(petri_net, computation="max", subspace=space)        
            
            
            # Go through every trap and see if it is already resolved.
            for trap in max_traps:
                #print(trap)
                trap_set = space_to_aeon_set(stg, trap)
                trap_candidates = trap_set.intersect(candidates)
                if not trap_candidates.is_empty():
                    attractors_recursive(trap, trap_candidates)
                    candidates = candidates.minus(trap_set)
                    all_traps = all_traps.union(trap_set)
                
        covered = reach_bwd(stg, all_traps, candidates)

        candidates = candidates.minus(covered)
        if not candidates.is_empty():
            # TODO: This is not correct, because it will look for sets of states that are
            # terminal within the candidates set. I.e. if there is a transition out from
            # candidates into a state that is not in candidates, this will not find that transition.
            symbolic_attractor_search(stg, candidates)
        
    attractors_recursive({}, stg.unit_colored_vertices())


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

    # Just a sanity check for .aeon models that can contain unknown update functions.
    if has_parameters(network):
        print("Network contains unknown parameters. This type of network is not supported yet.")
        sys.exit()

    # Disable static analysis in AEON (some randomly generated or 
    # pre-processed files might fail this check).
    network = remove_static_constraints(network)
    print(network)

    while True:
        attractors(network)
