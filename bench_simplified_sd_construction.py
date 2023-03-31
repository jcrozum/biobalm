from biodivine_aeon import BooleanNetwork
from nfvsmotifs.SuccessionDiagram import SuccessionDiagram
from nfvsmotifs.expansion_algorithms import simplified_bfs_expansion_2, simplified_bfs_expansion, simplified_dfs_expansion
import sys
import nfvsmotifs

# Print progress and succession diagram size.
nfvsmotifs.SuccessionDiagram.DEBUG = True
nfvsmotifs.expansion_algorithms.DEBUG = True

NODE_LIMIT = 1_000_000
DEPTH_LIMIT = 10_000

# This is unfortunately necessary for PyEDA Boolean expression parser (for now).
sys.setrecursionlimit(150000)

bn = BooleanNetwork.from_file(sys.argv[1])
bn = bn.infer_regulatory_graph()

# Compute the succession diagram.
sd = SuccessionDiagram(bn)
#expanded = simplified_bfs_expansion(sd, sd.root(), depth_limit=DEPTH_LIMIT, node_limit=NODE_LIMIT)
expanded = simplified_dfs_expansion(sd)

print(f"Succession diagram size:", sd.count_nodes())
print(f"Stub nodes:", sd.count_stubs())
print(f"Minimal traps:", len(sd.find_all_minimal_nodes()))