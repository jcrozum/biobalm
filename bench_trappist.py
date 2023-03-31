from biodivine_aeon import BooleanNetwork
from nfvsmotifs.trappist_core import trappist_async
import sys

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