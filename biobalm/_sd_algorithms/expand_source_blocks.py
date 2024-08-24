from __future__ import annotations

import itertools as it
from typing import TYPE_CHECKING, cast

from biobalm.types import BooleanSpace
from biobalm.interaction_graph_utils import source_nodes
from biobalm.space_utils import is_subspace, intersect

if TYPE_CHECKING:
    from biobalm.succession_diagram import SuccessionDiagram


def expand_source_blocks(
    sd: SuccessionDiagram,
    check_maa: bool = True,
    size_limit: int | None = None,
    optimize_source_nodes: bool = True,
    optimize_independent_blocks: bool = True,
) -> bool:
    """
    Base correctness assumptions:

     - Expanding two minimal blocks is always independent. Every minimal trap space is contained in one of the maximal trap
       spaces of a single block. Furthermore, expanding a minimal block cannot intrfere with the expansion of the other
       stable motifs.
     - Any path in the "inner" succession diagram based on a single block can be reproduced in the larger succession diagram,
       meaning absence of attractors in the inner diagram proves absence of attractors in the outer diagram.
    """

    if sd.config["debug"]:
        print(f"Start SD expansion using block decomposition on {sd.network}.")

    root = sd.root()

    current_level: set[int] = set([root])
    next_level: set[int] = set()

    bfs_depth = 0

    while len(current_level) > 0:
        bfs_depth += 1
        if sd.config["debug"]:
            total_fixed_vars = 0
            for node in current_level:
                total_fixed_vars += len(sd.node_data(node)["space"])
            print(
                f" > Start block expansion of a BFS level {bfs_depth} with {len(current_level)} node(s). Average number of fixed variables: {total_fixed_vars / len(current_level)}/{sd.network.variable_count()}"
            )

        for node in sorted(current_level):  # Sorted for determinism
            if sd.node_data(node)["expanded"]:
                # We re-discovered a previously expanded node.
                continue

            # Only continue if the succession diagram isn't too large.
            if (size_limit is not None) and (len(sd) >= size_limit):
                # Size limit reached.
                return False

            node_bn = sd.node_percolated_network(node, compute=True)
            node_space = sd.node_data(node)["space"]

            # Check for source nodes in the percolated network. We can "fast-forward" these the
            # same way we would skip them in the root node.
            sources = source_nodes(node_bn)
            if len(sources) != 0 and optimize_source_nodes:
                if sd.config["debug"]:
                    print(
                        f" > Found {len(sources)} source nodes in node {node}. Fast-forwarding instead of expansion."
                    )

                # Check if size limits are exceeded.
                expected_size = len(sd) + (2 ** len(sources))
                if expected_size > sd.config["max_motifs_per_node"]:
                    raise RuntimeError(
                        f"Exceeded the maximum amount of stable motifs per node ({sd.config['max_motifs_per_node']}; see `SuccessionDiagramConfiguration.max_motifs_per_node`)."
                    )
                elif size_limit is not None and expected_size > size_limit:
                    # Cannot expand. Size limit would be exceeded.
                    return False
                else:
                    if sd.config["debug"]:
                        print(
                            f" > Fast-forwarding {len(sources)} source nodes into {2 ** len(sources)} SD nodes."
                        )

                bin_values_iter = it.product(range(2), repeat=len(sources))
                for bin_values in bin_values_iter:
                    valuation = cast(BooleanSpace, dict(zip(sources, bin_values)))
                    sub_space = node_space | valuation

                    next_level.add(sd._ensure_node(node, sub_space))  # type: ignore

                # This makes the node artificially "expanded". Also, there
                # can be no attractors here because we are just fixing the source nodes.
                sd.node_data(node)["expanded"] = True
                sd.node_data(node)["attractor_seeds"] = []
                sd.node_data(node)["attractor_sets"] = []

                # This node is done, more work will be done at the next level.
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

            if len(minimal_blocks) > 1:

                if sd.config["debug"]:
                    print(
                        f" > [{node}] Try to fast-forward {len(minimal_blocks)} independent minimal blocks."
                    )

                all_min_traps: list[list[BooleanSpace]] = []
                sds = [
                    sd.component_subdiagram(list(x[0]), node) for x in minimal_blocks
                ]
                for i, s in enumerate(sds):
                    if (
                        s.network.variable_count() <= 3
                        or s.network.variable_count() > 100
                    ):
                        min_traps = [
                            sd.node_data(x)["space"] for x in minimal_blocks[i][1]
                        ]
                    else:
                        s.expand_minimal_spaces()
                        min_traps = [
                            (s.node_data(x)["space"] | node_space)
                            for x in s.minimal_trap_spaces()
                        ]
                        if sd.config["debug"]:
                            print(
                                f" > [{node}] Expanded {len(minimal_blocks[i][1])} block to {len(min_traps)} min traps."
                            )
                    all_min_traps.append(min_traps)

                created: set[int] = set()
                for combination in it.product(*all_min_traps, repeat=1):
                    combined_space: BooleanSpace = {}
                    for x in combination:
                        intersection = intersect(combined_space, x)
                        assert intersection is not None
                        combined_space = intersection

                    combined_node = sd._ensure_node(node, combined_space)  # type: ignore
                    created.add(combined_node)
                    # for x in combination:
                    sd._ensure_edge(node, combined_node, combined_space)  # type: ignore

                # created: set[int] = set()
                # block_sd_nodes = [ x[1] for x in minimal_blocks ]
                # for combination in it.product(*block_sd_nodes, repeat=1):
                #     combined_space: BooleanSpace = {}
                #     for x in combination:
                #         intersection = intersect(combined_space, sd.node_data(x)["space"])
                #         assert intersection is not None
                #         combined_space = intersection

                #     combined_node = sd._ensure_node(combination[0], combined_space) # type: ignore
                #     created.add(combined_node)
                #     for x in combination:
                #         sd._ensure_edge(x, combined_node, combined_space) # type: ignore

                next_level = next_level | created

                if sd.config["debug"]:
                    print(f" > [{node}] Fast-forwarded to nodes {created}.")

                continue

            # if len(minimal_blocks) > 1:
            #     # We have multiple independent minimal blocks. We can try to fast-forward the smallest one
            #     # to reduce the succession diagram size.

            #     block, block_nodes = minimal_blocks[0]
            #     block_sd = sd.component_subdiagram(list(block), node)

            #     if sd.config["debug"]:
            #         print(
            #             f" > [{node}] Try to fast-forward minimal block with {len(block)} variables."
            #         )

            #     assert block_sd.expand_block(find_motif_avoidant_attractors=True, optimize_source_nodes=True)

            #     if block_sd.node_data(block_sd.root())["attractor_seeds"] == []:
            #         sd.node_data(node)["attractor_seeds"] = []
            #         sd.node_data(node)["attractor_sets"] = []

            #     is_clean = True
            #     for inner_node in block_sd.expanded_ids():
            #         if block_sd.node_is_minimal(inner_node):
            #             continue

            #         inner_data = block_sd.node_data(inner_node)
            #         if inner_data["attractor_seeds"] != []:
            #             is_clean = False

            #     block_minimal_spaces = [ (block_sd.node_data(x)["space"] | node_space) for x in block_sd.minimal_trap_spaces() ]

            #     if sd.config["debug"]:
            #         print(
            #             f" > [{node}] Skipped {len(block_sd) - 1 - len(block_nodes)} nodes by fast-forwarding into {len(block_minimal_spaces)} minimal trap spaces."
            #         )

            #     for block_node in block_nodes:
            #         block_node_data = sd.node_data(block_node)
            #         if block_node_data["space"] in block_minimal_spaces:
            #             # This space is already in here.
            #             continue

            #         for successor in block_minimal_spaces:
            #             # Add links to all the *appropriate* minimal trap spaces.
            #             if is_subspace(successor, block_node_data["space"]):
            #                 sd._ensure_node(block_node, successor) # type: ignore

            #         block_node_data = sd.node_data(block_node)
            #         block_node_data["expanded"] = True

            #         if is_clean:
            #             block_node_data["attractor_seeds"] = []
            #             block_node_data["attractor_sets"] = []

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

                    # Instead of expanding the inner succession diagram "normally", we can
                    # copy the successors/stable motifs that we already know. This doesn't
                    # make too much of a difference unless the maximal trap space computation
                    # is complex, but that is sometimes true, especially for large networks.
                    for succ_id in block_nodes:
                        succ_motif = sd.edge_stable_motif(node, succ_id, reduced=True)
                        block_sd._ensure_node(block_sd.root(), succ_motif)  # type: ignore
                    block_sd.node_data(block_sd.root())["expanded"] = True

                    # We could also consider using `seeds` instead of `candidates` here. Ultimately, this
                    # matters very rarely. The reasoning for why we use `candidates` is that we can (almost)
                    # always guarantee that the expansion is going to finish, even if some nodes do not
                    # have their MAAs eliminated. As such, one can try to use other techniques to disprove
                    # MAAs in the problematic nodes while using the nice properties of the expansion to
                    # still disprove MAAs in the remaining nodes. If we used `seeds`, the expansion could
                    # just get stuck on this node and the "partial" results wouldn't be usable.
                    try:
                        block_sd_candidates = block_sd.node_attractor_candidates(
                            block_sd.root(), compute=True
                        )
                        is_clean = len(block_sd_candidates) == 0
                    except:
                        is_clean = False

                    if is_clean:
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
                                f"[{node}] > Found MAA candidates in a block or could not find candidates at all. Delaying expansion."
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
