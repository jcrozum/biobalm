from biodivine_aeon import BooleanNetwork, SymbolicAsyncGraph, find_attractors
from nfvsmotifs.SuccessionDiagram import SuccessionDiagram, simplying_dfs_expansion
from nfvsmotifs.sd_expansion import simplified_bfs_expansion
import sys
import nfvsmotifs

nfvsmotifs.SuccessionDiagram.DEBUG = True
nfvsmotifs.motif_avoidant.DEBUG = True
nfvsmotifs.sd_expansion.DEBUG = True
nfvsmotifs.sd_expansion.SYMBOLIC_PRUNING = True

NODE_LIMIT = 10_000

# This is unfortunately necessary for PyEDA Boolean expression parser (for now).
sys.setrecursionlimit(150000)

bn = BooleanNetwork.from_file(sys.argv[1])
bn = bn.infer_regulatory_graph()
stg = SymbolicAsyncGraph(bn)

# Compute the succession diagram.
sd = SuccessionDiagram(bn)
#simplying_dfs_expansion(sd)
actually_expanded = simplified_bfs_expansion(sd, sd.root())
#expanded = sd.expand_node(sd.root(), depth_limit=1000, node_limit=NODE_LIMIT)

#print("Expanded: ", expanded)
print("Nodes: {}", sd.G.number_of_nodes())
print("Expanded: {}", actually_expanded)
print("Stub nodes: ", sd.count_stubs())

# SD must be fully expanded, otherwise we may miss some results.
#assert expanded < NODE_LIMIT


# Compute attractors in diagram nodes.
# TODO: There will probably be a method that does this in one "go".
#nfvs_attractors = []
#for i in range(sd.G.number_of_nodes()):
#    attr = sd.expand_attractors(i)
#    for a in attr:
        # Just a simple sanity check.
#        assert len(a) == bn.num_vars()
#    if len(attr) > 0:
#        nfvs_attractors += attr

#print(f"Succession diagram attractors computed ({len(nfvs_attractors)} attractor[s]).")
#print("Start symbolic attractor computation...")

# Compute symbolic attractors using AEON.
#symbolic_attractors = find_attractors(stg)

# Check that every "seed" returned by SuccessionDiagram appears in
# some symbolic attractor, and that every symbolic attractor contains
# at most one such "seed" state.
#for seed in nfvs_attractors:
#    symbolic_seed = stg.fix_subspace({ k:bool(v) for k,v in seed.items() })
#    found = None
    
    # The "seed" state must have a symbolic attractor (and that
    # attractor mustn't have been removed yet).
#    for i in range(len(symbolic_attractors)):
#        if symbolic_seed.is_subset(symbolic_attractors[i]):
#            found = i
#    assert found is not None

#    symbolic_attractors.pop(found)

# All symbolic attractors must be covered by some seed at this point.
#assert len(symbolic_attractors) == 0

#print("All checks passed.")