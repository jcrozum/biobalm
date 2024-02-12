from __future__ import annotations

from typing import TYPE_CHECKING

from balm.symbolic_utils import state_list_to_bdd

if TYPE_CHECKING:
    from balm.SuccessionDiagram import SuccessionDiagram

import balm
import balm.SuccessionDiagram
from balm.motif_avoidant import detect_motif_avoidant_attractors, make_retained_set
from balm.trappist_core import compute_fixed_point_reduced_STG
from balm.types import BooleanSpace


def compute_attractor_seeds(
    sd: SuccessionDiagram,
    node_id: int,
) -> list[BooleanSpace]:
    """
    Compute the list of vertices such that each attractor within the subspace of
    the given `node_id` is covered by exactly one vertex.

    If the node is a stub, the result covers the whole subspace. If the node is
    expanded, the result only covers the "immediate" subspace without the
    subspaces of the child nodes.
    """

    if balm.SuccessionDiagram.DEBUG:
        print(f"[{node_id}] Start computing attractor seeds.")

    node_space = sd.node_space(node_id)

    if len(node_space) == sd.network.variable_count():
        # This node is a fixed-point.
        return [node_space]

    # Compute the list of child spaces if the node is expanded. Otherwise
    # "pretend" that there are no children.
    child_spaces = []
    if sd.node_is_expanded(node_id):
        child_spaces = [sd.node_space(s) for s in sd.node_successors(node_id)]

    # Fix everything in the NFVS to zero, as long as
    # it isn't already fixed by our `node_space`.
    #
    # We add the whole node space to the retain set because we know
    # the space is a trap and this will remove the corresponding unnecessary
    # Petri net transitions.
    retained_set = make_retained_set(
        sd.symbolic, sd.node_nfvs(node_id), node_space, child_spaces
    )

    if len(retained_set) == sd.network.variable_count() and len(child_spaces) == 0:
        # There is only a single attractor remaining here,
        # and its "seed" is the retained set.
        return [retained_set]

    terminal_restriction_space = state_list_to_bdd(
        sd.symbolic.symbolic_context(), child_spaces
    ).l_not()

    candidate_seeds = compute_fixed_point_reduced_STG(
        sd.petri_net,
        retained_set,
        ensure_subspace=node_space,
        avoid_subspaces=child_spaces,
    )

    if balm.SuccessionDiagram.DEBUG:
        print(f"[{node_id}] Found {len(candidate_seeds)} seed candidates.")

    if len(candidate_seeds) == 1 and len(child_spaces) == 0:
        # If this is a (non-strict) minimal trap and there is only one seed,
        # the seed must be valid.
        return candidate_seeds
    else:
        attractors = detect_motif_avoidant_attractors(
            sd.symbolic,
            sd.petri_net,
            candidate_seeds,
            terminal_restriction_space,
            max_iterations=1000,
            is_in_an_mts=len(child_spaces) == 0,
        )

        return attractors
