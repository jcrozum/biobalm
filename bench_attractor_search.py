from biodivine_aeon import BooleanNetwork
from nfvsmotifs.SuccessionDiagram import SuccessionDiagram
import sys
import nfvsmotifs

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
expanded = sd.expand_node(sd.root(), depth_limit=DEPTH_LIMIT, node_limit=NODE_LIMIT)

print("Succession diagram size:", expanded)

# SD must be fully expanded, otherwise we may miss some results.
assert expanded < NODE_LIMIT

# Compute attractors in succession diagram nodes.
nfvs_attractors = []
for i in range(sd.G.number_of_nodes()):
    attr = sd.expand_attractors(i)
    for a in attr:
        # Just a simple sanity check.
        assert len(a) == bn.num_vars()
    if len(attr) > 0:
        nfvs_attractors += attr

print("Attractors:", len(nfvs_attractors))