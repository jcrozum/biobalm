from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from biobalm.succession_diagram import SuccessionDiagram


def expand_bfs(
    sd: SuccessionDiagram,
    node_id: int | None = None,
    bfs_level_limit: int | None = None,
    size_limit: int | None = None,
) -> bool:
    """
    See `SuccessionDiagram.expand_bfs` for documentation.
    """

    if node_id is None:
        node_id = sd.root()

    seen: set[int] = set()
    seen.add(node_id)

    level_id = 0
    current_level = [node_id]
    next_level: list[int] = []

    while len(current_level) > 0:
        for node in current_level:
            # Check if the size limit has been exceeded already.
            if (size_limit is not None) and (len(sd) >= size_limit):
                # Size limit reached.
                return False

            # Compute successors if necessary.
            successors = sd.node_successors(node, compute=True)
            # Sort successors to avoid non-determinism.
            successors = sorted(successors)

            # Add successors to the next level and to the seen set.
            for s in successors:
                if s not in seen:
                    seen.add(s)
                    next_level.append(s)

        # The level is explored. Check if this exceeds the level limit.
        if (bfs_level_limit is not None) and (level_id >= bfs_level_limit):
            # Level limit reached.
            return False

        # If not, "move on" to the next level.
        level_id += 1
        current_level = next_level
        next_level = []

    return True
