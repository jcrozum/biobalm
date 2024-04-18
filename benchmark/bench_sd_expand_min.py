from biodivine_aeon import BooleanNetwork
from biobalm import SuccessionDiagram
import sys
import os
import pickle
import biobalm

# Print progress and succession diagram size.
biobalm.succession_diagram.DEBUG = True

NODE_LIMIT = 1_000_000
DEPTH_LIMIT = 10_000

bn = BooleanNetwork.from_file(sys.argv[1])
bn = bn.infer_valid_graph()

# Compute the succession diagram.
sd = SuccessionDiagram(bn)
sd.expand_minimal_spaces()

model_name = os.path.basename(sys.argv[1])
sd_name = os.path.splitext(model_name)[0] + ".pickle"
with open(f"./_sd_expand_min/{sd_name}", "wb") as handle:
    pickle.dump(sd, handle)

print(f"Succession diagram size:", len(sd))
print(f"Minimal traps:", len(sd.minimal_trap_spaces()))

print("size, expanded, minimal")
print(f"{len(sd)},{len(list(sd.expanded_ids()))},{len(sd.minimal_trap_spaces())}")
