import os
import sys
import time

from biodivine_aeon import BooleanNetwork # type:ignore
from networkx import DiGraph # type:ignore
from nfvsmotifs.interaction_graph_utils import infer_signed_interaction_graph, feedback_vertex_set, independent_cycles, find_minimum_NFVS

def main():
    model_dir = os.getcwd() + "/" + sys.argv[1]
    models = os.listdir(model_dir)
    models = sorted(models)
    
    for model in models:
            if not model.endswith(".bnet"):
                # Just in case there are some other files there.
                continue
            path = model_dir + "/" + model
            bn = BooleanNetwork.from_file(path)

            print(f"{model} ({len(bn.variables())} nodes)")

            start = time.process_time()
            nfvs_aeon = feedback_vertex_set(bn.infer_regulatory_graph().graph(), parity='negative')
            print(f"AEON: {len(nfvs_aeon)} ({time.process_time() - start:.2f})")
            
            start = time.process_time()
            nfvs_mtsNFVS = find_minimum_NFVS(bn)
            print(f"mtsNFVS: {len(nfvs_mtsNFVS)} ({time.process_time() - start:.2f})")

            print("=========")

if __name__ == "__main__":
    main()



