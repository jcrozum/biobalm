from biodivine_aeon import BooleanNetwork
from balm.trappist_core import trappist_async
import sys

# A simple script to benchmark minimal trap space detection using trappist.
# The output should be comparable to the output of `bench_sd_construction.py`
# in terms of minimal trap spaces.
#
# There is only one command line argument: a path to the network file.

sys.setrecursionlimit(150000)

bn = BooleanNetwork.from_file(sys.argv[1])
bn = bn.infer_regulatory_graph()

COUNT = 0

def count_all(x):
    global COUNT
    COUNT += 1
    return True

trappist_async(bn, on_solution=count_all, problem="min")

print("Minimal traps:", COUNT)

print("traps")
print(f"{COUNT}")