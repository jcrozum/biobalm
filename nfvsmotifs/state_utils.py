from __future__ import annotations
"""
    Operations on states
"""

from pyeda.boolalg import boolfunc # type:ignore
from pyeda.boolalg.bdd import bddvar, expr2bdd, BinaryDecisionDiagram # type:ignore
from pyeda.boolalg.expr import expr # type:ignore

from typing import List, Dict # type: ignore

def state_2_bdd(state: dict[str, int]) -> BinaryDecisionDiagram:
    state_bdd = 1

    for node in state:
        value = state[node]
        node_bdd = bddvar(node)

        if value == 1:
            state_bdd = state_bdd & node_bdd

        if value == 0:
            state_bdd = state_bdd & ~node_bdd

    return state_bdd


def list_state_2_bdd(states: list[dict[str, int]]) -> BinaryDecisionDiagram:
    result_bdd = 0

    for state in states:
        result_bdd = result_bdd | state_2_bdd(state)

    return result_bdd


def eval_function(f: BinaryDecisionDiagram, state: dict[str, int]) -> int:
    if f.is_zero():
        return 0

    state_bdd_dict = {}
    for node in state:
        node_bdd = bddvar(node)
        state_bdd_dict[node_bdd] = state[node]

    f_strict = f.restrict(state_bdd_dict)

    if f_strict.is_one():
        return 1
    else:
        return 0


def is_member_bdd(state: dict[str, int], f: BinaryDecisionDiagram) -> bool:
    if f.is_zero():
        return False

    state_bdd_dict = {}
    for node in state:
        node_bdd = bddvar(node)
        state_bdd_dict[node_bdd] = state[node]

    f_strict = f.restrict(state_bdd_dict)

    return f_strict.is_one()
