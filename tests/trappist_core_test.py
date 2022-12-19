from biodivine_aeon import BooleanNetwork, SymbolicAsyncGraph, FixedPoints
from nfvsmotifs.trappist_core import trappist
from nfvsmotifs.aeon_utils import remove_static_constraints
from nfvsmotifs.space_utils import is_trap_space
import os
import sys

# TODO: Right now, this is necessary to correctly parse all models.
# In the future, we should ideally use a parser that does not have this problem.
sys.setrecursionlimit(100_000)

#def test_trappist_spaces():
    # Verify that all minimum and maximum trap spaces are indeed trap spaces.
    #for model in os.listdir("./bbm-bnet-inputs-true"):
        #if not model.endswith(".bnet"):
            # Just in case there are some other files there.
        #    continue
        #bn = BooleanNetwork.from_file(f"./bbm-bnet-inputs-true/{model}")
        #if bn.num_vars() > 100:
            # We should be able to complete all inputs, but it would
            # take quite a long time. This seems like it can still finish 
            # in a few seconds.
        #    continue
        
        #for trap in trappist(bn, problem="min"):
        #    assert is_trap_space(bn, trap), f"Failed on {model}: {trap} is not a trap space."
        
        #for trap in trappist(bn, problem="max"):
        #    assert is_trap_space(bn, trap), f"Failed on {model}: {trap} is not a trap space."
    

def test_trappist_fixed_points():
    # Verify that the fixed-points of the test models are the same
    # as when computing using BDDs.
    for model in os.listdir("./bbm-bnet-inputs-true"):
        if not model.endswith(".bnet"):
            # Just in case there are some other files there.
            continue
        bn = BooleanNetwork.from_file(f"./bbm-bnet-inputs-true/{model}")
        if bn.num_vars() > 100:
            # We should be able to complete all inputs, but it would
            # take quite a long time. This seems like it can still finish 
            # in a few seconds.
            continue

        bn = remove_static_constraints(bn)
        stg = SymbolicAsyncGraph(bn)
        symbolic_fixed_points = FixedPoints.symbolic(stg)
        trappist_fixed_points = trappist(bn, problem="fix")
        for fixed_point in trappist_fixed_points:
            # Convert trappist result to a symbolic singleton set.
            vertex = stg.fix_subspace({ x: bool(int(fixed_point[x])) for x in fixed_point })
            # Check that the fixed-point has been found, and remove it.
            assert vertex.is_subset(symbolic_fixed_points)
            symbolic_fixed_points = symbolic_fixed_points.minus(vertex)
        # In the end, all fixed-points must have been found.
        assert symbolic_fixed_points.is_empty()        