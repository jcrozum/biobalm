"""
Optimized for finding the size of the fully expanded SD.
"""

import sys
sys.path.append("..")

import balm

from biodivine_aeon import BooleanNetwork
from balm.SuccessionDiagram import SuccessionDiagram


LOG_LOCATION = "SD_analysis_size_" + sys.argv[1].split("/")[-2] + ".csv"

# Print progress and succession diagram size.
balm.SuccessionDiagram.DEBUG = True # type: ignore

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


print(f"N,size,depth,att,maa,min_trap")
print(f"{bn.num_vars()},{len(sd)},{sd.depth()},,,{len(sd.minimal_trap_spaces())}")

log = open(LOG_LOCATION, "a")
log.write(sys.argv[1].split("/")[-1] + ",")  # model name
log.write(str(bn.num_vars()) + ",")  # network size
log.write(str(len(sd)) + ",")  # SD size
log.write(str(sd.depth()) + ",")  # SD depth
log.write(str(len(sd.minimal_trap_spaces())) + "\n")  # number of minimal trapspaces
