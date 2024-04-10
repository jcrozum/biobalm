from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from biobalm.succession_diagram import SuccessionDiagram

from biobalm.space_utils import is_subspace
from biobalm.trappist_core import trappist


def expand_minimal_spaces(sd: SuccessionDiagram, size_limit: int | None = None) -> bool:
    """
    See `SuccessionDiagram.expand_minimal_spaces` for documentation.
    """

    minimal_traps = trappist(sd.petri_net, problem="min")

    root = sd.root()

    seen = set([root])

    stack: list[tuple[int, list[int] | None]] = [(root, None)]

    while len(stack) > 0:
        (node, successors) = stack.pop()
        if successors is None:
            # Only allow successor computation if size limit hasn't been exceeded.
            if (size_limit is not None) and (len(sd) >= size_limit):
                # Size limit reached.
                return False

            successors = sd.node_successors(node, compute=True)
            successors = sorted(successors, reverse=True)  # For determinism!
            # (reversed because we explore the list from the back)

        node_space = sd.node_data(node)["space"]

        # Remove all immediate successors that are already visited or those who
        # do not cover any new minimal trap space.
        while len(successors) > 0:
            if successors[-1] in seen:
                successors.pop()
                continue
            if len([s for s in minimal_traps if is_subspace(s, node_space)]) == 0:
                successors.pop()
                continue
            break

        # This node is done and we don't have to push anything onto the stack.
        # In particular, it means that every minimal trap space which is a subset
        # of this node is already in the succession diagram.
        if len(successors) == 0:
            if sd.node_is_minimal(node):
                minimal_traps.remove(sd.node_data(node)["space"])
            continue

        # At this point, we know that `s` is not visited and it contains
        # at least one minimal trap space that does not appear in the
        # succession diagram yet.

        s = successors.pop()
        seen.add(s)
        # Push the node back with the remaining successors.
        stack.append((node, successors))
        # Push the successor onto the stack.
        stack.append((s, None))

    assert len(minimal_traps) == 0
    return True
