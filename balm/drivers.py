from __future__ import annotations

from biodivine_aeon import AsynchronousGraph

from balm.space_utils import percolate_space_strict
from balm.types import BooleanSpace


def find_single_node_LDOIs(stg: AsynchronousGraph) -> dict[tuple[str, int], BooleanSpace]:
    """
    finds LDOIs of every single node state
    TODO: take an initial set of LDOIs (e.g., of the original system) as an argument for speed-up
    """
    LDOIs: dict[tuple[str, int], BooleanSpace] = {}
    for var in stg.network_variable_names():
        fn_bdd = stg.mk_update_function(var)
        if fn_bdd.is_true() or fn_bdd.is_false():
            # Skip constant nodes.
            continue
        
        LDOIs[(var, 0)] = percolate_space_strict(stg, {var: 0})
        LDOIs[(var, 1)] = percolate_space_strict(stg, {var: 1})

    return LDOIs


def find_single_drivers(
    target_subspace: BooleanSpace,
    stg: AsynchronousGraph,
    LDOIs: dict[tuple[str, int], BooleanSpace] | None = None,
) -> set[tuple[str, int]]:
    """
    find all the single node drivers for a given target_subspace,
    usually (but not necessarily) a maximal trapspace (stablemotif)
    """
    if LDOIs is None:
        LDOIs = find_single_node_LDOIs(stg)

    drivers: set[tuple[str, int]] = set()
    for fix, LDOI in LDOIs.items():
        if target_subspace.items() <= (LDOI.items() | {fix}):
            drivers.add(fix)

    return drivers
