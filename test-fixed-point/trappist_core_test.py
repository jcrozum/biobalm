from biodivine_aeon import BooleanNetwork, SymbolicAsyncGraph, FixedPoints # type: ignore
from nfvsmotifs.trappist_core import trappist, compute_fixed_point_reduced_STG
from nfvsmotifs.aeon_utils import remove_static_constraints
from nfvsmotifs.space_utils import is_syntactic_trap_space
from nfvsmotifs.petri_net_translation import network_to_petrinet
import os
import sys

def test_network_fixed_point_reduced_STG():
    # Verify that the function for computing fixed points of the reduced STG

    bn = BooleanNetwork.from_bnet("""
        x1, (x1 & x2) | (!x1 & !x2)
        x2, (x1 & x2) | (!x1 & !x2)
    """)

    petri_net = network_to_petrinet(bn)
    nodes = []
    avoid_subspace_1 = {"x1" : "1", "x2" : "1"}
    avoid_subspace_2 = {}
    avoid_subspace_3 = {"x2" : "1"}

    ensure_subspace_1 = {}
    ensure_subspace_2 = {"x1" : "0", "x2" : "0"}

    for var in bn.variables():
        name = bn.get_variable_name(var)
        nodes.append(name)


    retained_set = {"x1" : "0", "x2" : "0"}
    candidate_set = compute_fixed_point_reduced_STG(petri_net, nodes, retained_set)
    assert len(candidate_set) == 2 # candidate_set = {00, 11}

    candidate_set = compute_fixed_point_reduced_STG(petri_net, nodes, retained_set, ensure_subspace=ensure_subspace_1, avoid_subspaces=[avoid_subspace_1])
    assert len(candidate_set) == 1 # candidate_set = {00}

    candidate_set = compute_fixed_point_reduced_STG(petri_net, nodes, retained_set, ensure_subspace=ensure_subspace_1, avoid_subspaces=[avoid_subspace_2])
    assert len(candidate_set) == 0 # candidate_set = empty

    candidate_set = compute_fixed_point_reduced_STG(petri_net, nodes, retained_set, ensure_subspace=ensure_subspace_2)
    assert len(candidate_set) == 1 # candidate_set = {00}

    retained_set = {"x1" : "1", "x2" : "1"}
    candidate_set = compute_fixed_point_reduced_STG(petri_net, nodes, retained_set)
    assert len(candidate_set) == 3 # candidate_set = {01, 10, 11}

    candidate_set = compute_fixed_point_reduced_STG(petri_net, nodes, retained_set, ensure_subspace=ensure_subspace_1, avoid_subspaces=[avoid_subspace_1])
    assert len(candidate_set) == 2 # candidate_set = {01, 10}

    candidate_set = compute_fixed_point_reduced_STG(petri_net, nodes, retained_set, ensure_subspace=ensure_subspace_1, avoid_subspaces=[avoid_subspace_2])
    assert len(candidate_set) == 0 # candidate_set = empty

    candidate_set = compute_fixed_point_reduced_STG(petri_net, nodes, retained_set, ensure_subspace=ensure_subspace_1, avoid_subspaces=[avoid_subspace_3])
    assert len(candidate_set) == 1 # candidate_set = {10}


    



