"""
    Utility operations on states. State (or a subspace) is represented as a dictionary of model
    variables mapping keys to 0/1 values.
"""

from __future__ import annotations
from balm.types import BooleanSpace

# TODO: 
#   For now, I deleted the code related to caching, because it is a bit harder to do when we need to also
#   manage the `SymbolicContext`. However, we should test that this does not impact performance.

def dnf_function_is_true(dnf: list[BooleanSpace], state: BooleanSpace) -> bool:
    """
    Returns `True` if the given DNF function evaluates to `1` for the given
    state (or space).
    """
    if len(dnf) == 0:
        return False

    for conjunction in dnf:
        if conjunction.items() <= state.items():
            return True
    return False


def remove_state_from_dnf(
    dnf: list[BooleanSpace], state: BooleanSpace
) -> list[BooleanSpace]:
    """
    Removes all conjunctions that are True in the state
    """
    modified_dnf: list[BooleanSpace] = []
    for conjunction in dnf:
        if conjunction.items() <= state.items():
            pass
        else:
            modified_dnf.append(conjunction)
    return modified_dnf
