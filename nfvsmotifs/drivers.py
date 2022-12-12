from __future__ import annotations

from biodivine_aeon import BooleanNetwork # type: ignore
from nfvsmotifs.space_utils import percolate_space

def find_single_node_LDOIs(bn: BooleanNetwork) -> dict[tuple[str, str], dict[str, str]]:
    """
    find LDOIs of single node fixes
    TODO: take LDOIs of the original system for speed-up
    """
    LDOIs = {}
    for var in bn.variables():
        name = bn.get_variable_name(var)
        function = bn.get_update_function(var)
        # TODO: other constant check
        if function == "true" or function == "false":
            continue
        for i in range(2):
            fix = (name, str(i))
            space = {name: str(i)}
            LDOIs[fix] = percolate_space(bn,space)[0]

    return LDOIs

def find_single_drivers(max_trap, bn: BooleanNetwork, LDOIs: dict = {}) -> list[tuple[str, str]]:
    """
    find all the single node drivers for a given maximal trapspace(stablemotif)
    """
    if LDOIs == {}:
        LDOIs = find_single_node_LDOIs(bn)

    drivers = []
    for fix in LDOIs:
        LDOI = LDOIs[fix]
        if max_trap.items() <= LDOI.items():
            drivers.append(fix)

    return drivers
