"""
Utility operations for creating and manipulating symbolic functions.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal, cast

if TYPE_CHECKING:
    from biodivine_aeon import Bdd, SymbolicContext

from biodivine_aeon import BddVariableSet, BddValuation

from biobalm.types import BooleanSpace


def state_to_bdd(
    bdd_context: SymbolicContext | BddVariableSet, state: BooleanSpace
) -> Bdd:
    """
    Convert a `BooleanSpace` to a BDD encoding the state singleton.

    Parameters
    ----------
    bdd_context : SymbolicContext | BddVariableSet
        The BDD context in which the state is represented. This is a
        a `biodivine_aeon.BddVariableSet` or `biodivine_aeon.SymbolicContext`
        (which is automatically converted to `biodivine_aeon.BddVariableSet`).
        The context object ensures compatibility between BDDs by
        maintaining a shared collection of variable names and their ordering.
    state : BooleanSpace
        The state to convert to a BDD.

    Returns
    -------
    Bdd
        The BDD encoding the given state.
    """
    bdd = (
        bdd_context
        if isinstance(bdd_context, BddVariableSet)
        else bdd_context.bdd_variable_set()
    )
    return bdd.mk_conjunctive_clause(state)


def state_list_to_bdd(
    bdd_context: SymbolicContext | BddVariableSet, states: list[BooleanSpace]
) -> Bdd:
    """
    Convert a list of state dictionaries to a BDD representation.

    The BDD will be satisfied in any state in the list and in no others.

    Parameters
    ----------
    bdd_context : SymbolicContext | BddVariableSet
        The BDD context in which the state is represented. This is a
        a `biodivine_aeon.BddVariableSet` or `biodivine_aeon.SymbolicContext`
        (which is automatically converted to `biodivine_aeon.BddVariableSet`).
        The context object ensures compatibility between BDDs by
        maintaining a shared collection of variable names and their ordering.
    states : list[BooleanSpace]
        The list of states to convert to a BDD.

    Returns
    -------
    Bdd
        The BDD encoding the given list of states.
    """
    bdd = (
        bdd_context
        if isinstance(bdd_context, BddVariableSet)
        else bdd_context.bdd_variable_set()
    )
    return bdd.mk_dnf(states)


def function_eval(f: Bdd, state: BooleanSpace) -> Literal[0, 1] | None:
    """
    Evaluate a BDD Boolean function in the given state.

    If the state is incomplete (i.e. it is a space), the function value can be
    undetermined. In such case, `None` is returned.

    Parameters
    ----------
    f : Bdd
        The function to evaluate.
    state : BooleanSpace
        The state in which the function is evaluated.

    Returns
    -------
    Literal[0, 1] | None
        The function value in the given state, or `None` if the function is
        undetermined.
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
    `True` if the BDD function evaluates to `1` for the given state (or space).

    Parameters
    ----------
    f : Bdd
        The function to evaluate.
    state : BooleanSpace
        The state in which the function is evaluated.

    Returns
    -------
    bool
        `True` if the function evaluates to `1` in the given state.
    """
    return function_eval(f, state) == 1


def valuation_to_state(ctx: SymbolicContext, valuation: BddValuation) -> BooleanSpace:
    """
    Extract network state from a `BddValuation` into a `BooleanSpace`.

    Parameters
    ----------
    f : SymbolicContext
        A context which maps between network variables and their symbolic counterparts.
    state : BddValuation
        A valuation of the network's symbolic encoding.

    Returns
    -------
    BooleanSpace
        A `BooleanSpace` encoding the state data from the given valuation.
    """

    result: BooleanSpace = {}
    for var, val in valuation.items():
        n_var = ctx.find_network_variable(var)
        if n_var is not None:
            result[ctx.get_network_variable_name(n_var)] = cast(Literal[0, 1], int(val))

    return result
