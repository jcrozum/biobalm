from biodivine_aeon import BooleanNetwork
from biobalm.petri_net_translation import network_to_petrinet
import biobalm.petri_net_translation
import sys

biobalm.petri_net_translation.DEBUG = True

bn = BooleanNetwork.from_file(sys.argv[1])
bn = bn.infer_valid_graph()

pn = network_to_petrinet(bn)

print("pn_nodes, pn_edges")
print(f"{len(pn.nodes())}, {len(pn.edges())}")
