from biodivine_aeon import BooleanNetwork, SymbolicAsyncGraph, find_attractors
import sys

bn = BooleanNetwork.from_file(sys.argv[1])
bn = bn.infer_regulatory_graph()

stg = SymbolicAsyncGraph(bn)

attr = find_attractors(stg)

print("Attractors:", len(attr))