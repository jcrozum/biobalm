import sys
sys.path.append(".")

from biodivine_aeon import BooleanNetwork
from nfvsmotifs.SuccessionDiagram import SuccessionDiagram
from nfvsmotifs._sd_algorithms.expand_source_SCCs import expand_source_SCCs
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
fully_expanded = expand_source_SCCs(sd)
assert fully_expanded

print(f"Succession diagram size:", len(sd))
print(f"Minimal traps:", len(sd.minimal_trap_spaces()))

print("size, expanded, minimal")
print(f"{len(sd)}, {len(list(sd.expanded_ids()))}, {len(sd.minimal_trap_spaces())}")