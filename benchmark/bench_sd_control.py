from biodivine_aeon import BooleanNetwork
from balm.SuccessionDiagram import SuccessionDiagram
from balm.control import succession_control
import sys
import time
import balm.SuccessionDiagram

# Print progress and succession diagram size.
balm.SuccessionDiagram.DEBUG = True

NODE_LIMIT = 1_000_000
DEPTH_LIMIT = 10_000

bn = BooleanNetwork.from_file(sys.argv[1])
bn = bn.infer_valid_graph()

# Compute the succession diagram.
sd = SuccessionDiagram(bn)
t0 = time.perf_counter()
fully_expanded = sd.expand_bfs(bfs_level_limit=DEPTH_LIMIT, size_limit=NODE_LIMIT)
t_expand = time.perf_counter() - t0
assert fully_expanded

# Compute control into the first minimal trap space of the diagram.
id = sd.minimal_trap_spaces()[0]
target = sd.node_space(id)

t0 = time.perf_counter()
interventions_nfvsmotifs2 = succession_control(bn, target, succession_diagram=sd)
t_control = time.perf_counter() - t0

print("nodes, interventions, control-runtime, expansion-runtime")
print(f"{len(sd)}, {len(interventions_nfvsmotifs2)}, {t_control}, {t_expand}")