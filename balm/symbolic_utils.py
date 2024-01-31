"""
    Utility operations for creating and manipulating symbolic functions.
"""

from __future__ import annotations
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Literal
    from biodivine_aeon import Bdd, SymbolicContext

from balm.types import BooleanSpace
from biodivine_aeon import BddVariableSet

def state_to_bdd(ctx: SymbolicContext | BddVariableSet, state: BooleanSpace) -> Bdd:
    """
    Convert a `BooleanSpace` to a BDD encoding the state singleton.
    """
    bdd = ctx if isinstance(ctx, BddVariableSet) else ctx.bdd_variable_set()
    return bdd.mk_conjunctive_clause(state)


def state_list_to_bdd(ctx: SymbolicContext | BddVariableSet, states: list[BooleanSpace]) -> Bdd:
    """
    Convert a list of state dictionaries to a BDD representation.
    """
    bdd = ctx if isinstance(ctx, BddVariableSet) else ctx.bdd_variable_set()
    return bdd.mk_dnf(states)


def function_eval(f: Bdd, state: BooleanSpace) -> Literal[0, 1] | None:
    """
    Evaluate a BDD function in the given state to an integer value. If the state is incomplete
    (i.e. it is a space), the function may not evaluate to an exact integer. In such case,
    `None` is returned.
    """
    if f.is_false():
        return 0
    if f.is_true():
        return 1    

    reduced_f = f.r_restrict(state)
    if reduced_f.is_true():
        return 1
    if reduced_f.is_false():
        return 0
    return None


def function_is_true(f: Bdd, state: BooleanSpace) -> bool:
    """
    Returns `True` if the given BDD function evaluates to `1` for the given
    state (or space).
    """
    return function_eval(f, state) == 1    
