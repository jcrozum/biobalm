from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from biobalm.succession_diagram import SuccessionDiagram

from biobalm.space_utils import intersect, is_subspace
from biobalm.types import BooleanSpace


def expand_to_target(
    sd: SuccessionDiagram, target: BooleanSpace, size_limit: int | None = None
):
    """
    See `SuccessionDiagram.exapnd_to_target` for documentation.
    """

    root = sd.root()
    seen = set([root])

    level_id = 0
    current_level = [root]
    next_level: list[int] = []

    while len(current_level) > 0:
        for node in current_level:
            node_space = sd.node_data(node)["space"]

            if intersect(node_space, target) is None:
                # If `node_space` does not intersect with `target`, it is not relevant
                # and we can safely keep it as unexpanded "stub".
                continue

            if is_subspace(node_space, target) and not node_space == target:
                # If `node_space` is a subspace of `target`, it is relevant but expanding
                # it will not add any new information, we can thus also keep it unexpanded.
                continue

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

        # If not, "move on" to the next level.
        level_id += 1
        current_level = next_level
        next_level = []

    return True
