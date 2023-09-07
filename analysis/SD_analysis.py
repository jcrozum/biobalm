import sys
sys.path.append("..")

from biodivine_aeon import BooleanNetwork
from nfvsmotifs.SuccessionDiagram import SuccessionDiagram
import nfvsmotifs

LOG_LOCATION = 'SD_analysis_'+sys.argv[1].split('/')[-2]+'.csv'

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
assert fully_expanded

attractor_count = 0
motif_avoidant_count = 0

for node in sd.node_ids():
    attr = sd.node_attractor_seeds(node, compute=True)
    attractor_count += len(attr)
    if not sd.node_is_minimal(node):
        motif_avoidant_count += len(attr)

print(f"network size, {bn.num_vars()}")
print(f"SD size, {len(sd)}")
print(f"SD depth, {sd.depth()} ")
print(f"attractors, {attractor_count}")
print(f"motif avoidant attactors, {motif_avoidant_count}")
print(f"minimal trapspaces, {len(sd.minimal_trap_spaces())}")

log = open(LOG_LOCATION, 'a')
log.write(sys.argv[1].split('/')[-1] + ',')             # model name
log.write(str(bn.num_vars()) + ',')                     # network size
log.write(str(len(sd)) + ',')                           # SD size
log.write(str(sd.depth()) + ',')                        # SD depth
log.write(str(attractor_count) + ',')                   # number of attractors
log.write(str(motif_avoidant_count) + ',')              # number of motif avoidant attractors
log.write(str(len(sd.minimal_trap_spaces())) + '\n')    # number of minimal trapspaces

