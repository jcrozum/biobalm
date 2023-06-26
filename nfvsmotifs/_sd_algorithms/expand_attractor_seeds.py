from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nfvsmotifs.SuccessionDiagram import SuccessionDiagram

from nfvsmotifs.trappist_core import compute_fixed_point_reduced_STG
from nfvsmotifs.motif_avoidant import make_retained_set
from nfvsmotifs.space_utils import intersect
import nfvsmotifs

def expand_attractor_seeds(sd: SuccessionDiagram):
    """
    See `SuccessionDiagram.expand_attractor_seeds` for documentation.
    """

    # First, expand the succession diagram such that all minimal trap spaces are found.
    # This reduces the amount of work performed in this algorithm, because for every attractor
    # in a minimal trap space, we already have the closest trap space, now we just need to
    # do the same for (potential) motif-avoidant attractors.
    sd.expand_minimal_spaces()

    if nfvsmotifs.SuccessionDiagram.DEBUG:
        print("Minimal trap space expansion finished. Proceeding to attractor expansion.")

    root = sd.root()
    seen = set([root])
    stack: list[tuple[int, list[int] | None]] = [(root, None)]

    while len(stack) > 0:
        (node, successors) = stack.pop()
        if successors is None:
            successors = sd.node_successors(node, compute=True)
            successors = sorted(successors, reverse=True) # For determinism!
            # (reversed because we explore the list from the back)
        
        node_space = sd.node_space(node)

        # Retrieve the stable motifs of children that are already expanded.
        expanded_children = [ x for x in sd.node_successors(node) if sd.node_is_expanded(x) ]
        expanded_motifs = [ sd.edge_stable_motif(node, child) for child in expanded_children ]

        # Now, we skip all successors that are either already seen, or that 
        # do not contain any candidate states for motif-avoidant attractors.
        while len(successors) > 0:
            if successors[-1] in seen:
                # The next node was already seen on stack. We can thus skip it and continue
                # to the next one.
                successors.pop()
                continue
            if sd.node_is_expanded(successors[-1]):
                # The next node to explore is expanded (by some previous procedure)
                # but not "seen" in this search yet. We need to visit this node
                # regardless of other conditions
                break
            # Now, we need to asses if the next successor has some candidate states which
            # are not covered by the already expanded children.

            successor_space = sd.node_space(successors[-1])
            retained_set = make_retained_set(sd.network, sd.nfvs, successor_space)

            avoid_or_none = [ intersect(successor_space, child) for child in expanded_motifs ]
            avoid = [ x for x in avoid_or_none if x is not None ]

            successor_seeds = compute_fixed_point_reduced_STG(
                sd.petri_net,
                retained_set,
                ensure_subspace=successor_space,
                avoid_subspaces=avoid,
                solution_limit=1,
            )

            if len(successor_seeds) == 0:
                # At this point, we know that this successor is not expanded and there are either
                # no candidate states in it, or all candidate states are already covered by some
                # other expanded successor.
                successors.pop()
                continue

            print(len(successor_seeds))

            if nfvsmotifs.SuccessionDiagram.DEBUG:
                print(f"[{node}] Found successor with new attractor candidate seeds. Expand node {successors[-1]}.")

            break

        if len(successors) == 0:
            # Everything is done for this `node` and we can continue to the next one.
            continue

        s = successors.pop()
        seen.add(s)
        # Push the node back with the remaining successors.
        stack.append((node, successors))
        # Push the successor onto the stack.
        stack.append((s, None))
            