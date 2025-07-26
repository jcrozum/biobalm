from biodivine_aeon import BooleanNetwork
from biobalm import SuccessionDiagram
from biobalm.types import BooleanSpace
from biobalm._sd_attractors.long_lived import compute_long_lived_phenotypes
import sys

# This script demonstrates a non-trivial attractor analysis of a Boolean network.
# It demonstrates possible configuration options for the attractor detection process.
# Alternatively, one can simply use `SuccessionDiagram.build()` and `SuccessionDiagram.summary()`
# to run attractor detection with default settings.

bn = BooleanNetwork.from_file(sys.argv[1])
print(f"Loaded network: {bn}")
print("Infering valid regulatory graph from update functions...")
bn = bn.infer_valid_graph()
print("Eliminating constant values...")
bn = bn.inline_constants(infer_constants=True, repair_graph=True)
print(f"Simplified network: {bn}")

# Prepare a config which will print progress.
config = SuccessionDiagram.default_config()
config["debug"] = True  # Print progress.
config["max_motifs_per_node"] = (
    1_000_000  # Maximum number of outgoing edges for each node.
)
config["attractor_candidates_limit"] = (
    100_000  # Maximum number of enumerated attractor candidates.
)

sd = SuccessionDiagram(bn, config)

# Compute full succession diagram:
sd.expand_dfs()

print(f"Succession diagram size:", len(sd))
print(f"Minimal traps:", len(sd.minimal_trap_spaces()))

result: dict[int, list[BooleanSpace]] = {}
all_long_lived: list[BooleanSpace] = []
for node in reversed(list(sd.expanded_ids())):
    long_lived = compute_long_lived_phenotypes(sd, node)
    node_data = sd.node_data(node)
    result[node] = long_lived
    for l in long_lived:
        if l not in all_long_lived:
            all_long_lived.append(l)


for node in reversed(list(sd.expanded_ids())):
    long_lived = result[node]
    node_data = sd.node_data(node)
    print(f"[ {node} | Fixed vars: {len(node_data['space'])} | {node_data['space']} ]")
    print(
        f"\t[ Long-lived: {len(long_lived)} | Fixed vars: {sorted([len(x) for x in long_lived])}]"
    )

print(f"Unique long-lived: {len(all_long_lived)}")
