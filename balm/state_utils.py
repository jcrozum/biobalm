from __future__ import annotations

from functools import cache

"""
    Utility operations on states. State (or a subspace) is represented as a dictionary of model
    variables mapping keys to 0/1 values.
"""

from typing import TYPE_CHECKING, Literal

from pyeda.boolalg.bdd import BDDONE, BDDZERO, bddvar  # type:ignore

if TYPE_CHECKING:
    from pyeda.boolalg.bdd import BDDVariable, BinaryDecisionDiagram

from balm.types import space_type


def _state_dict_to_bdd_valuation(state: space_type) -> dict[BDDVariable, int]:
    """
    Convert state variables in a dictionary to their BDD counterparts.
    """
    return {bddvar_cache(x): y for (x, y) in state.items()}


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


def state_to_bdd(state: space_type, usecache: bool = True) -> BinaryDecisionDiagram:
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


def state_list_to_bdd(states: list[space_type]) -> BinaryDecisionDiagram:
    """
    Convert a list of state dictionaries to a BDD representation.
    """
    result_bdd: BinaryDecisionDiagram = BDDZERO
    for state in states:
        result_bdd = result_bdd | state_to_bdd(state)

    return result_bdd


def function_restrict(
    f: BinaryDecisionDiagram, state: space_type
) -> BinaryDecisionDiagram:
    """
    Restrict the validity of the given BDD function to valuations which
    are compatible with the given state variable dictionary.
    """
    bdd_state = _state_dict_to_bdd_valuation(state)
    return f.restrict(bdd_state)


def function_eval(f: BinaryDecisionDiagram, state: space_type) -> Literal[0, 1] | None:
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


def function_is_true(f: BinaryDecisionDiagram, state: space_type) -> bool:
    """
    Returns `True` if the given BDD function evaluates to `1` for the given
    state (or space).
    """
    if f.is_zero():
        return False

    return function_restrict(f, state).is_one()


def dnf_function_is_true(dnf: list[space_type], state: space_type) -> bool:
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


def remove_state_from_dnf(dnf: list[space_type], state: space_type) -> list[space_type]:
    """
    Removes all conjunctions that are True in the state
    """
    modified_dnf: list[space_type] = []
    for conjunction in dnf:
        if conjunction.items() <= state.items():
            pass
        else:
            modified_dnf.append(conjunction)
    return modified_dnf
