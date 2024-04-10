from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from biobalm.succession_diagram import SuccessionDiagram


def expand_dfs(
    sd: SuccessionDiagram,
    node_id: int | None = None,
    dfs_stack_limit: int | None = None,
    size_limit: int | None = None,
) -> bool:
    """
    See `SuccessionDiagram.expand_dfs` for documentation.
    """

    if node_id is None:
        node_id = sd.root()

    seen: set[int] = set()
    seen.add(node_id)

    stack: list[tuple[int, list[int] | None]] = [(node_id, None)]

    result_is_complete = True

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

        # Remove all immediate successors that are already visited.
        while len(successors) > 0 and successors[-1] in seen:
            successors.pop()

        # This node is done and we don't have to push anything onto the stack.
        if len(successors) == 0:
            continue

        if (dfs_stack_limit is not None) and (len(stack) >= dfs_stack_limit):
            # We cannot push any successor nodes because it would exceed
            # the stack limit. As such, we can just continue with the next
            # item on the stack. however, we must remember that we skipped
            # some nodes and the result is thus incomplete.
            result_is_complete = False
            continue

        s = successors.pop()
        seen.add(s)
        # Push the node back with the remaining successors.
        stack.append((node, successors))
        # Push the successor onto the stack.
        stack.append((s, None))

    return result_is_complete
