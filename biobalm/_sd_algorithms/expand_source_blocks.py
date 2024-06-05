from __future__ import annotations

import itertools as it
from typing import TYPE_CHECKING, cast

from biobalm.space_utils import percolate_network
from biobalm.types import BooleanSpace
from biobalm.interaction_graph_utils import source_nodes

if TYPE_CHECKING:
    from biobalm.succession_diagram import SuccessionDiagram


def expand_source_blocks(
    sd: SuccessionDiagram,
    check_maa: bool = True,
) -> bool:
    """
    Base correctness assumptions:

     - Expanding two minimal blocks is always independent.
     -
    """

    if sd.config["debug"]:
        print(f"Start SD expansion using block decomposition on {sd.network}.")

    root = sd.root()

    current_level: set[int] = set([root])
    next_level: set[int] = set()

    # This already accounts for constant percolation.
    node_space = sd.node_data(root)["space"]

    # find source nodes
    perc_bn = percolate_network(sd.network, node_space)
    sources = source_nodes(perc_bn)

    if sd.config["debug"]:
        print(f" > Computed source/input variable(s): {sources}")

    # get source nodes combinations and expand root node
    if len(sources) != 0:
        # If there are too many source nodes, this can generate an absurdly large SD.
        # This would be a problem even without the SCC expansion, but we can just
        # stop the whole thing faster because we know how many nodes it generates.
        if 2 ** len(sources) > sd.config["max_motifs_per_node"]:
            raise RuntimeError(
                f"Exceeded the maximum amount of stable motifs per node ({sd.config['max_motifs_per_node']}; see `SuccessionDiagramConfiguration.max_motifs_per_node`)."
            )
        else:
            if sd.config["debug"]:
                print(
                    f" > Expanding {len(sources)} source node into {2 ** len(sources)} SD nodes."
                )

        bin_values_iter = it.product(range(2), repeat=len(sources))
        for bin_values in bin_values_iter:
            valuation = cast(BooleanSpace, dict(zip(sources, bin_values)))
            sub_space = node_space | valuation

            next_level.add(sd._ensure_node(root, sub_space))  # type: ignore

        # This makes the root artificially "expanded". Also, there
        # can be no attractors here because we are just fixing the source nodes.
        sd.node_data(root)["expanded"] = True
        sd.node_data(root)["attractor_seeds"] = []
        sd.node_data(root)["attractor_sets"] = []
        current_level = next_level
        next_level = set()

    bfs_depth = 0

    while len(current_level) > 0:
        bfs_depth += 1
        if sd.config["debug"]:
            print(
                f" > Start block expansion of a BFS level {bfs_depth} with {len(current_level)} node(s)."
            )

        for node in sorted(current_level):  # Sorted for determinism
            if sd.node_data(node)["expanded"]:
                # We re-discovered a previously expanded node.
                continue

            # Compute successors as if this was a normal expansion procedure.
            successors = sd.node_successors(node, compute=True)
            # Sort successors to avoid non-determinism.
            successors = sorted(successors)

            if len(successors) == 0:
                # This is a minimal trap space.
                continue

            if len(successors) == 1 and not check_maa:
                # This space is not minimal, but there is no "choice" to make here.
                # We can just continue with that single choice.
                #
                # However, if we are checking for MAAs, we continue as normal, because
                # we need to compute most of the non-trivial results anyway.

                next_level.add(successors[0])
                continue

            node_bn = sd.node_percolated_network(node, compute=True)

            # Split successors into "blocks" based on the regulatory component
            # of the variables fixed by the stable motif.

            # Maps a "block" (bwd-closed set of variables) to a list of node IDs (successor nodes).
            blocks: list[tuple[set[str], list[int]]] = []
            for s in successors:
                motif = sd.edge_stable_motif(node, s, reduced=True)
                motif_block = node_bn.backward_reachable(list(motif.keys()))
                motif_block_names = {node_bn.get_variable_name(v) for v in motif_block}

                found = False
                for block, nodes in blocks:
                    if block == motif_block_names:
                        found = True
                        nodes.append(s)
                        break
                if not found:
                    blocks.append((motif_block_names, [s]))

            if sd.config["debug"]:
                print(
                    f" > [{node}] Computed blocks: {[(len(k), len(v)) for (k, v) in blocks]}"
                )

            # Now remove all non-minimal blocks (minimal in terms of set inclusion).
            # The reason why we are removing those is that they are not independent
            # on the minimal blocks.

            minimal_blocks: list[tuple[set[str], list[int]]] = []
            if len(blocks) > 1:
                for block, nodes in blocks:
                    is_minimal = True
                    for b2, _ in blocks:
                        if b2 < block:
                            is_minimal = False
                            break
                    if is_minimal:
                        minimal_blocks.append((block, nodes))
            else:
                minimal_blocks = blocks

            # Sort minimal blocks by the number of successor nodes (we want to pick
            # the one that leads leads to the minimal SD expansion).
            minimal_blocks = sorted(minimal_blocks, key=lambda x: len(x[1]))

            if sd.config["debug"]:
                print(
                    f" > [{node}] Minimal blocks: {[(len(k), len(v)) for (k, v) in minimal_blocks]}"
                )

            if not check_maa:
                # We will expand all nodes that are in the smallest block.
                to_expand = minimal_blocks[0][1]

                if sd.config["debug"]:
                    print(f" > [{node}] Final block ({len(to_expand)}): {to_expand}")

                next_level = next_level | set(to_expand)
            else:
                # Here, we want to find the smallest block without any MAAs and choose it.
                # If such block does not exist, we expand the whole node, because the MAAs
                # can be either in this node, or in any of the child nodes. If a clean block
                # is found, we know that it is safe to expand it without "missing" any MAAs.
                clean_block_found = False
                for block, block_nodes in minimal_blocks:
                    block_sd = sd.component_subdiagram(list(block), node)

                    # The succession diagram "restricted" to the considered block should have
                    # the same (restricted) successor nodes.
                    assert len(
                        block_sd.node_successors(block_sd.root(), compute=True)
                    ) == len(block_nodes)

                    # We could also consider using `seeds` instead of `candidates` here. Ultimately, this
                    # matters very rarely. The reasoning for why we use `candidates` is that we can (almost)
                    # always guarantee that the expansion is going to finish, even if some nodes do not
                    # have their MAAs eliminated. As such, one can try to use other techniques to disprove
                    # MAAs in the problematic nodes while using the nice properties of the expansion to
                    # still disprove MAAs in the remaining nodes. If we used `seeds`, the expansion could
                    # just get stuck on this node and the "partial" results wouldn't be usable.
                    block_sd_candidates = block_sd.node_attractor_candidates(
                        block_sd.root(), compute=True
                    )
                    if len(block_sd_candidates) == 0:
                        if sd.config["debug"]:
                            print(
                                f" > [{node}] Found clean block with no MAAs ({len(block_nodes)}): {block_nodes}"
                            )
                        clean_block_found = True
                        next_level = next_level | set(block_nodes)
                        sd.node_data(node)["attractor_seeds"] = []
                        sd.node_data(node)["attractor_sets"] = []
                        break
                    else:
                        if sd.config["debug"]:
                            print(
                                f"[{node}] > Found {len(block_sd_candidates)} MAA cnadidates in a block. Delaying expansion."
                            )
                if not clean_block_found:
                    # If all blocks have MAAs, we expand all successors.
                    if sd.config["debug"]:
                        print(
                            f" > [{node}] No clean block found. Expanding all {len(successors)} successors."
                        )
                    next_level = next_level | set(successors)

        current_level = next_level
        next_level = set()

    if sd.config["debug"]:
        print(f" > Block expansion terminated with {len(sd)} node(s).")

    return True
