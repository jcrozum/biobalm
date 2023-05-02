from __future__ import annotations
from functools import cache
"""
    Utility operations on states. State (or a subspace) is represented as a dictionary of model 
    variables mapping keys to 0/1 values.
"""

from pyeda.boolalg.bdd import bddvar, BDDONE, BDDZERO # type:ignore

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from pyeda.boolalg.bdd import BinaryDecisionDiagram, BDDVariable

def _state_dict_to_bdd_valuation(state: dict[str, int]) -> dict[BDDVariable, int]:
    """
        Convert state variables in a dictionary to their BDD counterparts.
    """
    return { bddvar_cache(x): y for (x,y) in state.items() }

@cache
def bddvar_cache(name: str) -> BDDVariable:
    return bddvar(name)

@cache
def state_to_bdd_cacheable(state: frozenset[tuple[str, int]]) -> BinaryDecisionDiagram:
    """
        Convert a state variables to a BDD encoding the state singleton.
    """
    state_bdd: BinaryDecisionDiagram = BDDONE

    for node, value in state:
        node_bdd = bddvar_cache(node)

        if value == 1:
            state_bdd = state_bdd & node_bdd

        if value == 0:
            state_bdd = state_bdd & ~node_bdd

    return state_bdd

def state_to_bdd(state: dict[str, int], usecache: bool = True) -> BinaryDecisionDiagram:
    """
        Convert a state variables to a BDD encoding the state singleton.
    """
    if usecache:
        return state_to_bdd_cacheable(frozenset(state.items()))
    state_bdd: BinaryDecisionDiagram = BDDONE

    for node, value in state.items():
        node_bdd = bddvar(node)

        if value == 1:
            state_bdd = state_bdd & node_bdd

        if value == 0:
            state_bdd = state_bdd & ~node_bdd

    return state_bdd


def state_list_to_bdd(states: list[dict[str, int]]) -> BinaryDecisionDiagram:
    """
        Convert a list of state dictionaries to a BDD representation.
    """
    result_bdd: BinaryDecisionDiagram = BDDZERO
    for state in states:
        result_bdd = result_bdd | state_to_bdd(state)

    return result_bdd

def function_restrict(f: BinaryDecisionDiagram, state: dict[str, int]) -> BinaryDecisionDiagram:
    """
        Restrict the validity of the given BDD function to valuations which 
        are compatible with the given state variable dictionary.
    """
    bdd_state = _state_dict_to_bdd_valuation(state)
    return f.restrict(bdd_state)

def function_eval(f: BinaryDecisionDiagram, state: dict[str, int]) -> int | None:
    """
        Evaluate a BDD function in the given state to an integer value. If the state is incomplete
        (i.e. it is a space), the function may not evaluate to an exact integer. In such case,
        `None` is returned.
    """
    if f.is_zero():
        return 0

    reduced_f = function_restrict(f, state)
    if reduced_f.is_one():
        return 1
    if reduced_f.is_zero():
        return 0
    return None

def function_is_true(f: BinaryDecisionDiagram, state: dict[str, int]) -> bool:
    """
        Returns `True` if the given BDD function evaluates to `1` for the given 
        state (or space). 
    """
    if f.is_zero():
        return False

    return function_restrict(f, state).is_one()

def dnf_function_is_true(dnf: list[dict[str, int]], state: dict[str, int]) -> bool:
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

def remove_state_from_dnf(dnf: list[dict[str, int]], state: dict[str, int]) -> list[dict[str, int]]:
    """
        Removes all conjunctions that are True in the state
    """
    modified_dnf: list[dict[str, int]] = []
    for conjunction in dnf:
        if conjunction.items() <= state.items():
            pass
        else:
            modified_dnf.append(conjunction)
    return modified_dnf