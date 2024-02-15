from biodivine_aeon import BooleanNetwork
from balm.trappist_core import trappist_async
from balm.types import BooleanSpace
import sys

# A simple script to benchmark minimal trap space detection using trappist.
# The output should be comparable to the output of `bench_sd_construction.py`
# in terms of minimal trap spaces.
#
# There is only one command line argument: a path to the network file.

bn = BooleanNetwork.from_file(sys.argv[1])
bn = bn.infer_valid_graph()

count = 0


def count_all(x: BooleanSpace):
    global count
    count += 1
    return True


trappist_async(bn, on_solution=count_all, problem="min")

print("Minimal traps:", count)

print("traps")
print(f"{count}")
