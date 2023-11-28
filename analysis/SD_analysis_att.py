"""
Optimized for finding the attractors
"""

import sys
sys.path.append("..")

import nfvsmotifs

from biodivine_aeon import BooleanNetwork
from nfvsmotifs.SuccessionDiagram import SuccessionDiagram
from nfvsmotifs._sd_algorithms.expand_source_SCCs import expand_source_SCCs


LOG_LOCATION = "SD_analysis_att_" + sys.argv[1].split("/")[-2] + ".csv"

# Print progress and succession diagram size.
nfvsmotifs.SuccessionDiagram.DEBUG = True # type: ignore

# This is unfortunately necessary for PyEDA Boolean expression parser (for now).
sys.setrecursionlimit(150000)


bn = BooleanNetwork.from_file(sys.argv[1])
bn = bn.infer_regulatory_graph()

# Compute the succession diagram.
sd = SuccessionDiagram(bn)
fully_expanded = expand_source_SCCs(sd)
assert fully_expanded

# Find the attractors
attractor_count = 0
motif_avoidant_count = 0

for node in sd.node_ids():
    attr = sd.node_attractor_seeds(node, compute=True)
    attractor_count += len(attr)
    if not sd.node_is_minimal(node):
        motif_avoidant_count += len(attr)


print(f"N,size,depth,att,maa,min_trap")
print(f"{bn.num_vars()},{len(sd)},{sd.depth()},{attractor_count},{motif_avoidant_count},{len(sd.minimal_trap_spaces())}")

log = open(LOG_LOCATION, "a")
log.write(sys.argv[1].split("/")[-1] + ",")  # model name
log.write(str(bn.num_vars()) + ",")  # network size
log.write(str(sd.depth()) + ",")  # SD depth
log.write(str(attractor_count) + ",")  # number of attractors
log.write(str(motif_avoidant_count) + ",")  # number of motif avoidant attractors
log.write(str(len(sd.minimal_trap_spaces())) + "\n")  # number of minimal trapspaces