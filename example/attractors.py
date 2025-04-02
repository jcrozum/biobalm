from biodivine_aeon import BooleanNetwork
from biobalm import SuccessionDiagram
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
config[
    "max_motifs_per_node"
] = 1_000_000  # Maximum number of outgoing edges for each node.
config[
    "attractor_candidates_limit"
] = 100_000  # Maximum number of enumerated attractor candidates.

sd = SuccessionDiagram(bn, config)

# Compute block-based expansion of the succession diagram.
fully_expanded = sd.expand_block(
    find_motif_avoidant_attractors=True,
    exact_attractor_detection=True,
    size_limit=1_000,  # Maximum size of the succession diagram.
)

# If the succession diagram is too large to create fully, mark
# unexpanded nodes as skipped and continue with the analysis.
# Attractors in skip nodes are still found, but the same attractor can
# sometimes be found multiple times in different nodes.
skipped = 0
if not fully_expanded:
    skipped = sd.skip_remaining()

print(f"Succession diagram size:", len(sd))
print(f"Minimal traps:", len(sd.minimal_trap_spaces()))
if skipped > 0:
    print(f"Skip nodes: {skipped}")

attractor_count = 0
motif_avoidant_count = 0

# The order of iteration is important, because it ensures we start
# from minimal (or rather last expanded) nodes, i.e. the info from these
# can be reused in the other methods.
for node in reversed(list(sd.expanded_ids())):
    attr = sd.node_attractor_seeds(node, compute=True, symbolic_fallback=True)
    print(
        f"Found {len(attr)} attractor(s) in node {node}. List of attractor seeds: {attr}"
    )
    attractor_count += len(attr)
    if len(attr) > 0 and not sd.node_is_minimal(node):
        print("Node is not minimal. Motif-avoidant attractor found.")
        motif_avoidant_count += len(attr)
