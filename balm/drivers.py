from __future__ import annotations

from typing import Literal, cast

from biodivine_aeon import BooleanNetwork

from balm.space_utils import percolate_space
from balm.types import BooleanSpace


def find_single_node_LDOIs(bn: BooleanNetwork) -> dict[tuple[str, int], BooleanSpace]:
    """
    finds LDOIs of every single node state
    TODO: take an initial set of LDOIs (e.g., of the original system) as an argument for speed-up
    """
    LDOIs: dict[tuple[str, int], BooleanSpace] = {}
    for var in bn.variables():
        name = bn.get_variable_name(var)
        function = bn.get_update_function(var)
        # TODO: other constant check
        if function == "true" or function == "false":
            continue
        for i in range(2):
            fix = (name, i)
            space: BooleanSpace = {name: cast(Literal[0, 1], i)}
            LDOIs[fix] = percolate_space(bn, space)

    return LDOIs


def find_single_drivers(
    target_subspace: BooleanSpace,
    bn: BooleanNetwork,
    LDOIs: dict[tuple[str, int], BooleanSpace] | None = None,
) -> set[tuple[str, int]]:
    """
    find all the single node drivers for a given target_subspace,
    usually (but not necessarily) a maximal trapspace (stablemotif)
    """
    if LDOIs is None:
        LDOIs = find_single_node_LDOIs(bn)

    drivers: set[tuple[str, int]] = set()
    for fix, LDOI in LDOIs.items():
        if target_subspace.items() <= (LDOI.items() | {fix}):
            drivers.add(fix)

    return drivers
