from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from biobalm.succession_diagram import SuccessionDiagram

import copy
from biobalm.space_utils import is_subspace
from biobalm.trappist_core import trappist
from biobalm.types import BooleanSpace


def expand_minimal_spaces(
    sd: SuccessionDiagram,
    node_id: int | None,
    size_limit: int | None = None,
    skip_remaining: bool = False,
) -> bool:
    """
    See `SuccessionDiagram.expand_minimal_spaces` for documentation.
    """

    if node_id is None:
        node_id = sd.root()

    pn = sd.node_percolated_petri_net(node_id, compute=True)
    node_space = sd.node_data(node_id)["space"]

    all_minimal_traps = trappist(network=pn, problem="min", ensure_subspace=node_space)
    all_minimal_traps = [(node_space | x) for x in all_minimal_traps]
    # We don't need to duplicate the actual trap spaces, just the list.
    minimal_traps = copy.copy(all_minimal_traps)

    if sd.config["debug"]:
        print(
            f"Start minimal trap SD expansion using {len(minimal_traps)} minimal traps."
        )

    seen = set([node_id])

    stack: list[tuple[int, list[int] | None]] = [(node_id, None)]

    def make_skip_node(
        sd: SuccessionDiagram, node_id: int, all_minimal_traps: list[BooleanSpace]
    ):
        node = sd.node_data(node_id)
        if node["expanded"]:
            return

        skip_edges = 0
        for m_trap in all_minimal_traps:
            if is_subspace(m_trap, sd.node_data(node_id)["space"]):
                m_id = sd._ensure_node(node_id, m_trap)  # type: ignore
                m_data = sd.node_data(m_id)
                m_data["expanded"] = True
                assert sd.node_is_minimal(m_id)
                skip_edges += 1

        node["expanded"] = True
        node["skipped"] = True

        if sd.config["debug"]:
            print(f"[{node_id}] Node skipped with {skip_edges} edges.")

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
                # Everything in seen is expanded, so no need to skip it.
                successors.pop()
                continue
            if len([s for s in minimal_traps if is_subspace(s, node_space)]) == 0:
                skipped = successors.pop()
                if skip_remaining:
                    make_skip_node(sd, skipped, all_minimal_traps)
                continue
            break

        # This node is done and we don't have to push anything onto the stack.
        # In particular, it means that every minimal trap space which is a subset
        # of this node is already in the succession diagram.
        if len(successors) == 0:
            if sd.node_is_minimal(node):
                minimal_traps.remove(sd.node_data(node)["space"])
                if sd.config["debug"]:
                    print(f"Remaining minimal traps: {len(minimal_traps)}.")
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

        if sd.config["debug"]:
            print(f"[{s}] Expanding...")

    assert len(minimal_traps) == 0
    return True
