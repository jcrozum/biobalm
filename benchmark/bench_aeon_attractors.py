from biodivine_aeon import BooleanNetwork, AsynchronousGraph, Attractors
import sys

# This is a very basic script for running attractor detection usin AEON.
# The results should be comparable to `bench_attractor_search`.
#
# The script only takes one argument: a path to the network file.

bn = BooleanNetwork.from_file(sys.argv[1])
bn = bn.infer_valid_graph()

stg = AsynchronousGraph(bn)

attr = Attractors.attractors(stg, stg.mk_unit_colored_vertices())

print("Attractors:", len(attr))

attr_states = stg.mk_empty_colored_vertices()

for a in attr:
    attr_states = attr_states.union(a)

print("attractors, attractor-states")
print(f"{len(attr)},{attr_states.cardinality()}")