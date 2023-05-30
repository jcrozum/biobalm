from biodivine_aeon import BooleanNetwork
from nfvsmotifs.SuccessionDiagram import SuccessionDiagram
import sys
import nfvsmotifs

# A simple script to benchmark attractor detection using succession diagrams.
# The script only has one argument: a path to the network file.
# Currently, the script performs full succession diagram expansion.

# Print progress and succession diagram size.
nfvsmotifs.SuccessionDiagram.DEBUG = True

NODE_LIMIT = 1_000_000
DEPTH_LIMIT = 10_000

# This is unfortunately necessary for PyEDA Boolean expression parser (for now).
sys.setrecursionlimit(150000)

bn = BooleanNetwork.from_file(sys.argv[1])
bn = bn.infer_regulatory_graph()

# Compute the succession diagram.
sd = SuccessionDiagram(bn)
fully_expanded = sd.expand_bfs(bfs_level_limit=DEPTH_LIMIT, size_limit=NODE_LIMIT)

print("Succession diagram size:", len(sd))

# SD must be fully expanded, otherwise we may miss some results.
assert fully_expanded

# Compute attractors in succession diagram nodes.
nfvs_attractors = []
for i in sd.node_ids():
    attr = sd.node_attractor_seeds(i, compute=True)
    for a in attr:
        # Just a simple sanity check.
        assert len(a) == bn.num_vars()
    if len(attr) > 0:
        nfvs_attractors += attr

print("Attractors:", len(nfvs_attractors))