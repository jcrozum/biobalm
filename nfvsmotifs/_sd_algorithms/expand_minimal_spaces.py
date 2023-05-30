from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from nfvsmotifs.SuccessionDiagram import SuccessionDiagram

    
from nfvsmotifs.space_utils import is_subspace
from nfvsmotifs.trappist_core import trappist

def expand_minimal_spaces(sd: SuccessionDiagram):
    """
    See `SuccessionDiagram.expand_minimal_spaces` for documentation.
    """

    minimal_traps = trappist(sd.petri_net, problem="min")

    root = sd.root()

    seen = set()
    seen.add(root)

    level_id = 0
    current_level = [root]
    next_level = []

    while len(current_level) > 0:
        remaining_traps = minimal_traps.copy()

        # First, eliminate all trap spaces that are covered by an already expanded node.
        for node in current_level:
            if not sd.node_is_expanded(node):
                continue

            node_space = sd.node_space(node)
            # Remove all traps that are subspaces of this expanded node.        
            remaining_traps = [ t for t in remaining_traps if not is_subspace(t, node_space)]

            # Also add all successors to the next level, we may need to explore these further.
            for s in sorted(sd.node_successors(node)):
                if s not in seen:
                    seen.add(s)
                    next_level.append(s)

        # Now, do the same thing again, but this time actually expand nodes.
        for node in current_level:
            if sd.node_is_expanded(node):
                continue    # Expanded nodes are already processed.

            node_space = sd.node_space(node)
            original = len(remaining_traps)
            remaining_traps = [ t for t in remaining_traps if not is_subspace(t, node_space)]
            updated = len(remaining_traps)

            if updated == original:
                # This node does not cover any new minimal trap spaces. It is thus safe to
                # keep it as an unexpanded stub.
                continue

            successors = sd.node_successors(node, compute=True)
            successors = sorted(successors)

            # Add successors to the next level and to the seen set.
            for s in successors:
                if s not in seen:
                    seen.add(s)
                    next_level.append(s)

        assert len(remaining_traps) == 0

        level_id += 1
        current_level = next_level
        next_level = []

    return True
