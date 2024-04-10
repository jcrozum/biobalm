"""
Some basic utility operations on spaces (partial assignments of BN variables).

Each space is represented as a dictionary with a subset of variable names as
keys and values `0`/`1` assigned to fixed variables. See also
:class:`BooleanSpace<biobalm.types.BooleanSpace>`.
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

from biobalm.symbolic_utils import function_eval

if TYPE_CHECKING:
    from biodivine_aeon import BooleanExpression

    from biobalm.types import BooleanSpace


def intersect(x: BooleanSpace, y: BooleanSpace) -> BooleanSpace | None:
    """
    Compute the intersection of two spaces.

    Returns the intersection as a new space, or `None` if the spaces don't
    intersect.

    Parameters
    ----------
    x : BooleanSpace
        The first space.
    y : BooleanSpace
        The second space.

    Returns
    -------
    BooleanSpace | None
        The intersection of `x` and `y`, or `None` if the spaces don't
        intersect.
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

    Parameters
    ----------
    x : BooleanSpace
        The first space.
    y : BooleanSpace
        The second space.

    Returns
    -------
    bool
        `True` if `x` is a subspace of `y`.
    """
    for var in y:
        if var not in x:
            return False
        if x[var] != y[var]:
            return False
    return True


def dnf_function_is_true(dnf: list[BooleanSpace], state: BooleanSpace) -> bool:
    """
    Checks if a DNF function evaluates to `1` for the given state (or space).

    The DNF function is represented as a list of spaces, such that each space represents
    exactly one DNF clause.

    Parameters
    ----------
    dnf : list[BooleanSpace]
        The DNF function to evaluate.
    state : BooleanSpace
        The state in which the function is evaluated.

    Returns
    -------
    bool
        `True` if the function evaluates to `1` in the given state.
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
    Removes all clauses (conjunctions) that are `True` in the given `state` from a DNF function.

    The result is a new list (i.e., it does not modify the original list).

    Parameters
    ----------
    dnf : list[BooleanSpace]
        The DNF function to modify.
    state : BooleanSpace
        The state to remove from the function.

    Returns
    -------
    list[BooleanSpace]
        The modified DNF function.
    """
    modified_dnf: list[BooleanSpace] = []
    for conjunction in dnf:
        if conjunction.items() <= state.items():
            pass
        else:
            modified_dnf.append(conjunction)
    return modified_dnf


def percolate_space_strict(
    network: AsynchronousGraph, space: BooleanSpace
) -> BooleanSpace:
    """
    Percolates a space through a Boolean network, disregarding constants.

    Returns the set of variables which become fixed as a result of fixing the
    variables from `space` within the given `AsynchronousGraph`.

    Note that the strict percolation process does not propagate constants that
    are already fixed within the `network`, only those that are specified in
    `space`. Also, the result only contains any *new* constants, not those that
    are already fixed in `space`.

    Parameters
    ----------
    network : AsynchronousGraph
        A symbolic representation of the Boolean network (from which state
        transitions can be generated) in which the percolation is performed. If
        needed, an `AsynchronousGraph` object can be constructed from a
        `BooleanNetwork` via `biodivine_aeon.AsynchronousGraph(bn)`.
    space : BooleanSpace
        The space to percolate.

    Returns
    -------
    BooleanSpace
        The percolated space.
    """

    result: BooleanSpace = {}
    restriction: BooleanSpace = copy(space)
    candidates = set(network.network_variable_names())

    # Ignore variables that are already fixed.
    for var in network.network_variable_names():
        fn_bdd = network.mk_update_function(var)
        if fn_bdd.is_true() or fn_bdd.is_false():
            candidates.remove(var)

    done = False
    while not done:
        done = True
        for var in copy(candidates):
            fn_bdd = network.mk_update_function(var)
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
    network: AsynchronousGraph,
    space: BooleanSpace,
) -> BooleanSpace:
    """
    Percolates a space through a Boolean network.

    Takes a symbolic `AsynchronousGraph` and a `BooleanSpace`. It then percolates
    any values that are effectively constant within the `network` assuming the variables
    from `space` are fixed accordingly.

    If the argument is a trap space, then the result is a subspace of the
    argument and is also a trap space.

    However, when the argument is a general space, the percolation can actually
    lead "outside" of the original space. In such case, the original fixed value
    is *not* modified and the conflict will remain in the resulting space.

    Parameters
    ----------
    network : AsynchronousGraph
        A symbolic representation of the Boolean network (from which state
        transitions can be generated) in which the percolation is performed. If
        needed, an `AsynchronousGraph` object can be constructed from a
        `BooleanNetwork` via `biodivine_aeon.AsynchronousGraph(bn)`.
    space : BooleanSpace
        The space to percolate.

    Returns
    -------
    BooleanSpace
        The percolated space.
    """

    percolated = Percolation.percolate_subspace(network, space)
    result: BooleanSpace = {}
    for var, value in percolated.items():
        var_name = network.get_network_variable_name(var)
        result[var_name] = cast(Literal[0, 1], int(value))
    return result


def percolation_conflicts(
    network: AsynchronousGraph,
    space: BooleanSpace,
    strict_percolation: bool = True,
) -> set[str]:
    """
    Find variables that conflict with the percolation of the given space.

    Returns a set of variables from `space` that are in conflict with the
    percolation of the given space (see `percolate_space`).

    Parameters
    ----------
    network : AsynchronousGraph
        A symbolic representation of the Boolean network (from which state
        transitions can be generated) in which the percolation is performed. If
        needed, an `AsynchronousGraph` object can be constructed from a
        `BooleanNetwork` via `biodivine_aeon.AsynchronousGraph(bn)`.
    space : BooleanSpace
        The space to percolate.
    strict_percolation : bool
        If `True` (the default), then the percolation is performed using
        :func:`percolate_space_strict`. Otherwise, it is performed using
        :func:`percolate_space`.

    Returns
    -------
    set[str]
        A set of variables from `space` that are in conflict with the
        percolation of the given space.
    """
    conflicts: set[str] = set()

    if strict_percolation:
        perc_space = percolate_space_strict(network, space)
    else:
        perc_space = percolate_space(network, space)

    for var, value in perc_space.items():
        fn_bdd = network.mk_update_function(var)
        fn_value = function_eval(fn_bdd, perc_space)
        if fn_value is not None and value != fn_value:
            conflicts.add(var)

    return conflicts


def percolate_network(
    bn: BooleanNetwork,
    space: BooleanSpace,
    symbolic_network: AsynchronousGraph | None = None,
    remove_constants: bool = False,
) -> BooleanNetwork:
    """
    Reduces a Boolean network by percolating a given space.

    Takes a `BooleanNetwork` and a `BooleanSpace`. It then produces a new network with
    update functions percolated based on the supplied space.

    The dynamics of the resulting network correspond to the dynamics of the
    network obtained by percolating the given space. If the space (or the
    percolated space) is a trap space, then the resulting dynamics are a
    subgraph of the original network's state transition graph. Otherwise, the
    dynamics correspond to the effects of an external intervention.

    The percolation process is based on BDD conversion. For this purpose, an optional
    `SymbolicContext` can be provided. If not given, a temporary `SymbolicContext` will
    be created instead. Note that this is necessary to resolve non-trivial tautologies or
    contradictions that can arise once the variables from `space` are fixed.

    Parameters
    ----------
    bn : BooleanNetwork
        The network to percolate.
    space : BooleanSpace
        The space to percolate.
    symbolic_network : AsynchronousGraph | None
        An optional symbolic representation to use to perform the percolation. If not
        given, a temporary one will be created from `bn`.
    remove_constants : bool
        If `True`, then the constants are removed from the resulting network. By
        default, `False`.

    Returns
    -------
    BooleanNetwork
        The percolated network.
    """

    if symbolic_network is None:
        symbolic_network = AsynchronousGraph(bn)
    var_set = symbolic_network.symbolic_context().bdd_variable_set()

    # Percolate the space first to ensure everything that can be fixed is fixed.
    space = percolate_space(symbolic_network, space)

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
            percolated = restrict_expression(
                update.as_expression(), space, symbolic_context=var_set
            )
            new_update = UpdateFunction(new_bn, percolated)
            new_bn.set_update_function(var, new_update)

    new_bn = new_bn.infer_valid_graph()
    if remove_constants:
        new_bn = new_bn.inline_constants(infer_constants=True, repair_graph=True)

    return new_bn


def restrict_expression(
    expression: BooleanExpression,
    space: BooleanSpace,
    symbolic_context: BddVariableSet | SymbolicContext | None = None,
) -> BooleanExpression:
    """
    Restricts a Boolean expression to a given space.

    Takes a `BooleanExpression` and a `BooleanSpace`. Returns a simplified `BooleanExpression`
    that is valid for exactly the same members of the given `space` as the original expression.
    The resulting expression does not depend on the variables which are fixed in the given `space`.

    The percolation process is based on BDD conversion. For this purpose, an optional `BddVariableSet`
    can be provided. If not given, a temporary `BddVariableSet` will be created instead. Note that
    this is necessary to resolve non-trivial tautologies/contradictions that can arise once the
    variables from `space` are fixed.

    Parameters
    ----------
    expression : BooleanExpression
        The expression to restrict.
    space : BooleanSpace
        The space to restrict to.
    symbolic_context : BddVariableSet | SymbolicContext | None
        An optional symbolic context to use to perform the percolation. If not given,
        a temporary one will be created.

        This is a `biodivine_aeon.BddVariableSet` or a `biodivine_aeon.SymbolicContext`
        (which is automatically converted to `biodivine_aeon.BddVariableSet`).
        The context object ensures compatibility between BDDs by
        maintaining a shared collection of variable names and their ordering.

    Returns
    -------
    BooleanExpression
        The restricted expression.
    """

    variables = expression.support_set()
    space = {k: v for k, v in space.items() if k in variables}

    if len(space) == 0:
        return expression

    if symbolic_context is None:
        symbolic_context = BddVariableSet(sorted(variables))
    if isinstance(symbolic_context, SymbolicContext):
        symbolic_context = symbolic_context.bdd_variable_set()

    bdd = symbolic_context.eval_expression(expression)
    bdd = bdd.r_restrict(space)
    return bdd.to_expression()


def expression_to_space_list(
    expression: BooleanExpression,
    symbolic_context: BddVariableSet | SymbolicContext | None = None,
) -> list[BooleanSpace]:
    """
    Convert a Boolean expression to a list of subspaces for which it is true.

    Equivalent to a disjunctive normal form. Convert a `BooleanExpression` to a
    list of subspaces whose union represents an equivalent set of the network
    states which satisfy the expression.

    Note that the spaces are not necessarily pair-wise disjoint. Also, the list
    is not necessarily minimal.

    The translation uses a DNF conversion based on BDDs. For this purpose, an optional
    `BddVariableSet` can be provided. If not given, a temporary `BddVariableSet` will be
    created instead.

    Parameters
    ----------
    expression : BooleanExpression
        The expression to convert.
    symbolic_context : BddVariableSet | SymbolicContext | None
        An optional symbolic context to use to perform the percolation. If not given,
        a temporary one will be created.

        This is a `biodivine_aeon.BddVariableSet` or a `biodivine_aeon.SymbolicContext`
        (which is automatically converted to `biodivine_aeon.BddVariableSet`).
        The context object ensures compatibility between BDDs by
        maintaining a shared collection of variable names and their ordering.

    Returns
    -------
    list[BooleanSpace]
        The list of subspaces on which the expression is true.
    """

    if symbolic_context is None:
        variables = sorted(expression.support_set())
        symbolic_context = BddVariableSet(variables)
    if isinstance(symbolic_context, SymbolicContext):
        symbolic_context = symbolic_context.bdd_variable_set()

    bdd = symbolic_context.eval_expression(expression)

    sub_spaces: list[BooleanSpace] = []
    for clause in bdd.clause_iterator():
        space: BooleanSpace = {}
        for var, value in clause.items():
            space[symbolic_context.get_variable_name(var)] = cast(
                Literal[0, 1], int(value)
            )
        sub_spaces.append(space)

    return sub_spaces


def space_unique_key(space: BooleanSpace, network: BooleanNetwork) -> int:
    """
    Provide a unique hash key for the provided space in a given network.

    Computes an integer which is a unique representation of the provided `space`
    (with respect to the given `network`).

    This integer key can be used instead of the original `space` in places where
    dictionaries are not allowed, such as a key within a larger dictionary, or
    a sorting key.

    Note that when used for sorting, this key essentially implements a particular
    form of lexicographic ordering on spaces. This is always a total ordering
    (there is no ambiguity).

    Parameters
    ----------
    space : BooleanSpace
        The space to encode.
    network : BooleanNetwork
        The network in which the space is defined.

    Returns
    -------
    int
        A unique key for the space.
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
