"""
    Some basic utility operations on spaces (partial assignments of BN variables).

    Each space is represented as a dictionary with a subset of variable names as
    keys and values `0`/`1` assigned to fixed variables.
"""
from __future__ import annotations

from copy import copy
from typing import TYPE_CHECKING, Literal, cast

from biodivine_aeon import (
    AsynchronousGraph,
    BddVariableSet,
    BooleanNetwork,
    Percolation,
    SymbolicContext,
    UpdateFunction,
)

from balm.symbolic_utils import function_eval

if TYPE_CHECKING:
    from biodivine_aeon import BooleanExpression

    from balm.types import BooleanSpace


def intersect(x: BooleanSpace, y: BooleanSpace) -> BooleanSpace | None:
    """
    Compute the space which is the intersection of two spaces, or `None` if the spaces
    don't intersect.
    """
    result: BooleanSpace = {}
    for k, v in x.items():
        result[k] = v
    for k, v in y.items():
        if k in result and result[k] != v:
            return None
        result[k] = v
    return result


def is_subspace(x: BooleanSpace, y: BooleanSpace) -> bool:
    """
    Checks if `x` is a subspace of `y`.
    """
    for var in y:
        if var not in x:
            return False
        if x[var] != y[var]:
            return False
    return True


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


def percolate_space_strict(stg: AsynchronousGraph, space: BooleanSpace) -> BooleanSpace:
    """
    Returns the set of variables which become fixed as a result of fixing the
    variables from `space` within the given `AsynchronousGraph`.

    Note that the strict percolation process does not propagate constants that
    are already fixed within the `stg`, only those that are specified in `space`.
    Also, the result only contains any *new* constants, not those that are already
    fixed in `space`.
    """

    result: BooleanSpace = {}
    restriction: BooleanSpace = copy(space)
    candidates = set(stg.network_variable_names())

    # Ignore variables that are already fixed.
    for var in stg.network_variable_names():
        fn_bdd = stg.mk_update_function(var)
        if fn_bdd.is_true() or fn_bdd.is_false():
            candidates.remove(var)

    done = False
    while not done:
        done = True
        for var in copy(candidates):
            fn_bdd = stg.mk_update_function(var)
            fn_value = function_eval(fn_bdd, restriction)
            if fn_value is not None:
                if var in restriction and restriction[var] != fn_value:
                    # There is a conflict. We don't want to output this,
                    # but we also don't want to change the value.
                    candidates.remove(var)
                else:
                    done = False
                    restriction[var] = fn_value
                    result[var] = fn_value
                    candidates.remove(var)

    return result


def percolate_space(
    stg: AsynchronousGraph,
    space: BooleanSpace,
) -> BooleanSpace:
    """
    Takes a symbolic `AsynchronousGraph` and a `BooleanSpace`. It then percolates
    any values that are effectively constant with the `stg` assuming the variables
    from `space` are fixed accordingly.

    If the argument is a trap space, then the result is a subspace of the
    argument and is also a trap space.

    However, when the argument is a general space, the percolation can actually
    lead "outside" of the original space. In such case, the original fixed value
    is *not* modified and the conflict will remain in the resulting space.
    """

    percolated = Percolation.percolate_subspace(stg, space)
    result: BooleanSpace = {}
    for var, value in percolated.items():
        var_name = stg.get_network_variable_name(var)
        result[var_name] = cast(Literal[0, 1], int(value))
    return result


def percolation_conflicts(
    stg: AsynchronousGraph,
    space: BooleanSpace,
    strict_percolation: bool = True,
) -> set[str]:
    """
    Returns a set of variables from `space` that are in conflict with the percolation of
    the given space (see `percolate_space`).
    """
    conflicts: set[str] = set()

    if strict_percolation:
        perc_space = percolate_space_strict(stg, space)
    else:
        perc_space = percolate_space(stg, space)

    for var, value in perc_space.items():
        fn_bdd = stg.mk_update_function(var)
        fn_value = function_eval(fn_bdd, perc_space)
        if fn_value is not None and value != fn_value:
            conflicts.add(var)

    return conflicts


def percolate_network(
    bn: BooleanNetwork, space: BooleanSpace, ctx: AsynchronousGraph | None = None
) -> BooleanNetwork:
    """
    Takes a `BooleanNetwork` and a `BooleanSpace`. It then produces a new network with
    update functions percolated based on the supplied space.

    There are two caveats to this operation:

        (1) If the given space is *not* a trap space, it is up to you to figure
        out what is the relationship between the dynamics of the original and
        the percolated network. For trap spaces, we know that every transition
        inside that trap space is preserved exactly.

    The percolation process is based on BDD conversion. For this purpose, an optional
    `SymbolicContext` can be provided. If not given, a temporary `SymbolicContext` will
    be created instead. Note that this is necessary to resolve non-trivial tautologies or
    contradictions that can arise once the variables from `space` are fixed.
    """

    if ctx is None:
        ctx = AsynchronousGraph(bn)
    var_set = ctx.symbolic_context().bdd_variable_set()

    # Percolate the space first to ensure everything that can be fixed is fixed.
    space = percolate_space(ctx, space)

    # Make a copy of the BN and copy the relevant functions.
    new_bn = copy(bn)

    for var in bn.variables():
        update = bn.get_update_function(var)
        if update is None:
            # This variable is a free input.
            assert len(bn.predecessors(var)) == 0
            name = bn.get_variable_name(var)
            if name in space:
                new_bn.set_update_function(
                    var, UpdateFunction.mk_const(new_bn, space[name])
                )
        else:
            percolated = percolate_expression(
                update.as_expression(), space, ctx=var_set
            )
            new_update = UpdateFunction(new_bn, percolated)
            new_bn.set_update_function(var, new_update)

    return new_bn.infer_valid_graph()


def percolate_expression(
    expression: BooleanExpression,
    space: BooleanSpace,
    ctx: BddVariableSet | SymbolicContext | None = None,
) -> BooleanExpression:
    """
    Takes a `BooleanExpression` and a `BooleanSpace`. Returns a simplified `BooleanExpression`
    that is valid for exactly the same members of the given `space` as the original expression.
    The resulting expression does not depend on the variables which are fixed in the given `space`.

    The percolation process is based on BDD conversion. For this purpose, an optional `BddVariableSet`
    can be provided. If not given, a temporary `BddVariableSet` will be created instead. Note that
    this is necessary to resolve non-trivial tautologies/contradictions that can arise once the
    variables from `space` are fixed.
    """

    variables = expression.support_set()
    space = {k: v for k, v in space.items() if k in variables}

    if len(space) == 0:
        return expression

    if ctx is None:
        ctx = BddVariableSet(sorted(variables))
    if isinstance(ctx, SymbolicContext):
        ctx = ctx.bdd_variable_set()

    bdd = ctx.eval_expression(expression)
    bdd = bdd.r_restrict(space)
    return bdd.to_expression()


def expression_to_space_list(
    expression: BooleanExpression, ctx: BddVariableSet | SymbolicContext | None = None
) -> list[BooleanSpace]:
    """
    Convert a `BooleanExpression` to a list of subspaces whose union represents
    an equivalent set of the network states which satisfy the expression.

    Note that the spaces are not necessarily pair-wise disjoint. Also, the list
    is not necessarily minimal.

    The translation uses a DNF conversion based on BDDs. For this purpose, an optional
    `BddVariableSet` can be provided. If not given, a temporary `BddVariableSet` will be
    created instead.
    """

    if ctx is None:
        variables = sorted(expression.support_set())
        ctx = BddVariableSet(variables)
    if isinstance(ctx, SymbolicContext):
        ctx = ctx.bdd_variable_set()

    bdd = ctx.eval_expression(expression)

    sub_spaces: list[BooleanSpace] = []
    for clause in bdd.clause_iterator():
        space: BooleanSpace = {}
        for var, value in clause.items():
            space[ctx.get_variable_name(var)] = cast(Literal[0, 1], int(value))
        sub_spaces.append(space)

    return sub_spaces


def space_unique_key(space: BooleanSpace, network: BooleanNetwork) -> int:
    """
    Computes an integer which is a unique representation of the provided `space`
    (with respect to the given `network`).

    This integer key can be used instead of the original `space` in places where
    dictionaries are not allowed, such as a key within a larger dictionary, or
    a sorting key.

    Note that when used for sorting, this key essentially implements a particular
    form of lexicographic ordering on spaces. This is always a total ordering
    (there is no ambiguity).
    """

    # Key is a binary encoding of the space dictionary. Since Python has
    # arbitrary-precision integers, this should work for any network and be
    # reasonably fast (we are not doing any copies or string manipulation).
    key: int = 0
    for k, v in space.items():
        var = network.find_variable(k)
        if var is None:
            raise IndexError(f"Unknown variable {var}.")
        # Each variable is encoded as two bits, so the total length
        # of the key is 2 * n and the offset of each variable is 2 * index.
        # 00 - unknown; 10 - zero; 11 - one
        key |= (v + 2) << (2 * int(var))
    return key
