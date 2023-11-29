from __future__ import annotations

from nfvsmotifs.state_utils import bddvar_cache, function_restrict

"""
    Some basic utility operations on spaces (partial assignments of BN variables).

    Each space is represented as a dictionary with a subset of variable names as
    keys and values `0`/`1` assigned to fixed variables.
"""

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from pyeda.boolalg.expr import Expression
    from pyeda.boolalg.bdd import BinaryDecisionDiagram

from biodivine_aeon import BooleanNetwork, RegulatoryGraph
from pyeda.boolalg.expr import Complement, Literal, Variable

from nfvsmotifs.pyeda_utils import (
    PYEDA_FALSE,
    PYEDA_TRUE,
    aeon_to_bdd,
    aeon_to_pyeda,
    pyeda_to_aeon,
    substitute_variables_in_expression,
)


def intersect(x: dict[str, int], y: dict[str, int]) -> dict[str, int] | None:
    """
    Compute the space which is the intersection of two spaces, or `None` if the spaces
    don't intersect.
    """
    result: dict[str, int] = {}
    for k, v in x.items():
        result[k] = v
    for k, v in y.items():
        if k in result and result[k] != v:
            return None
        result[k] = v
    return result


def is_subspace(x: dict[str, int], y: dict[str, int]) -> bool:
    """
    Checks if `x` is a subspace of `y`.
    """
    for var in y:
        if var not in x:
            return False
        if x[var] != y[var]:
            return False
    return True


def is_syntactic_trap_space(bn: BooleanNetwork, space: dict[str, int]) -> bool:
    """
    Uses percolation to check if the given `space` is a trap space in the given `BooleanNetwork`.

    Note that this does not perform any sophisticated "semantic" analysis of the update functions.
    For example, if the update function contains a contradiction/tautology, this method will not
    be able to take this into account. However, aside from such "degenerate" cases, this should
    still work in typical practical scenarios.

    If you need a guaranteed test, you can try `SymbolicAsyncGraph::is_trap_set` instead.
    """
    for var in bn.variables():
        var_name = bn.get_variable_name(var)

        if var_name in space:
            expression = aeon_to_pyeda(bn.get_update_function(var))
            expression = percolate_pyeda_expression(expression, space)
            if str(space[var_name]) != str(expression):
                print(
                    space[var_name], str(expression), bn.get_update_function(var), space
                )
                return False
    return True


def percolate_space(
    network: BooleanNetwork,
    space: dict[str, int],
    strict_percolation: bool = True,
) -> dict[str, int]:
    """
    Takes a Boolean network and a space (partial assignment of `0`/`1` to the
    network variables). It then percolates the values in the given `space` to
    the remaining network variables based on the update functions of the given
    `network`.

    If the argument is a trap space, then the result is a subspace of the
    argument and is also a trap space.

    However, when the argument is a general space, the percolation can actually
    lead "outside" of the original space. In such case, the original fixed value
    is *not* modified and the conflict will remain in the resulting space.

    If `strict_percolation` is used, only variables that become fixed as a
    result of fixing the space are considered (e.g., nodes with constant update
    functions are not propagated). Furthermore, the variables in `space` are
    only returned if the value of their update funciton becomes fixed as a
    result of percolating the fixed node values specified by `space`.
    """

    if strict_percolation:
        result: dict[str, int] = {}
    else:
        result = {var: space[var] for var in space}

    fixed = {var: space[var] for var in space}
    fixed_bddvars = {bddvar_cache(k): v for k, v in fixed.items()}
    bdds: dict[str, BinaryDecisionDiagram] = {}
    bdd_inputs = {}
    var_name_dict = {network.get_variable_name(var): var for var in network.variables()}
    deletion_list = list(result)

    done = False
    while not done:
        for k in deletion_list:
            del var_name_dict[k]
        deletion_list: list[str] = []
        done = True
        for var_name, var in var_name_dict.items():
            if var_name in result:
                continue
            if var_name not in bdds:
                bdds[var_name] = aeon_to_bdd(network.get_update_function(var))

            if bdds[var_name].is_one():
                if strict_percolation:
                    continue
                r = 1
            elif bdds[var_name].is_zero():
                if strict_percolation:
                    continue
                r = 0
            else:
                if var_name not in bdd_inputs:
                    bdd_inputs[var_name] = set(bdds[var_name].inputs)  # type: ignore
                point = {x: fixed_bddvars[x] for x in bdd_inputs[var_name] if x in fixed_bddvars}  # type: ignore
                if point:
                    bdds[var_name] = bdds[var_name].restrict(point)  # type: ignore
                    bdd_inputs[var_name] -= set(point)  # type: ignore

                if bdds[var_name].is_one():
                    r = 1
                elif bdds[var_name].is_zero():
                    r = 0
                else:
                    r = -1
                    continue

            if var_name not in fixed:
                fixed[var_name] = r
                fixed_bddvars[bddvar_cache(var_name)] = r
                result[var_name] = r
                deletion_list.append(var_name)
                done = False
            elif fixed[var_name] == r and var_name not in result:
                result[var_name] = r

    return result


def percolation_conflicts(
    network: BooleanNetwork,
    space: dict[str, int],
    strict_percolation: bool = True,
) -> set[str]:
    """
    Returns a set of variable names that are in conflict with the percolation of
    the given space (see `percolate_space`).
    """
    conflicts: set[str] = set()

    perc_space = percolate_space(network, space, strict_percolation)

    for var, value in space.items():
        function_restricted = function_restrict(
            aeon_to_bdd(network.get_update_function(var)), perc_space
        )

        if function_restricted.is_one() and value == 0:
            conflicts.add(var)
        elif function_restricted.is_zero() and value == 1:
            conflicts.add(var)

    return conflicts


def percolate_network(bn: BooleanNetwork, space: dict[str, int]) -> BooleanNetwork:
    """
    Takes an AEON.py Boolean network and a space (partial assignment of
    network variables to `0`/`1`). It then produces a new network with
    update functions percolated based on the supplied space.
    There are two caveats to this operation:

        (1) If the given space is *not* a trap space, it is up to you to figure
        out what the relationship between the original and the resulting dynamics
        is. For trap spaces, we know that everything inside that trap space
        is preserved. If the space is not a trap space, you have now cut away
        all outgoing transitions.
        (2) The underlying regulatory graph of the new network retains all
        regulations of the original network, but all integrity constraints
        (essentiality, monotonicity) are removed, because they most likely
        no longer hold in the new network.
    """
    # Make a copy of the original regulatory network, but without integrity constraints.
    old_rg = bn.graph()
    new_rg = RegulatoryGraph([bn.get_variable_name(var) for var in bn.variables()])

    for reg in old_rg.regulations():
        reg["observable"] = False
        if "monotonicity" in reg:
            del reg["monotonicity"]
        new_rg.add_regulation(reg)

    # Copy the Boolean network, but with simplified expressions.
    new_bn = BooleanNetwork(new_rg)
    for var in bn.variables():
        name = bn.get_variable_name(var)
        new_expr = None
        if name in space:
            # If the value is fixed, just use it as a value directly.
            if space[name] == 1:
                new_expr = PYEDA_TRUE
            elif space[name] == 0:
                new_expr = PYEDA_FALSE
            else:
                raise ValueError(f"{space[name]=} is not a valid variable value")
        else:
            # If the value is not fixed, use a simplified expression.
            expression = aeon_to_pyeda(bn.get_update_function(var))
            new_expr = percolate_pyeda_expression(expression, space)

        new_bn.set_update_function(var, pyeda_to_aeon(new_expr))

    return new_bn


def percolate_pyeda_expression(
    expression: Expression, space: dict[str, int]
) -> Expression:
    """
    Takes a PyEDA expression and a subspace (dictionary assigning `1`/`0` to
    a subset of variables). Returns a simplified expression that is valid
    for exactly the same members of the given `space` as the original expression.
    The resulting expression does not depend on the variables which are fixed
    in the given `space`.
    """
    substitution = {x: (PYEDA_TRUE if space[x] == 1 else PYEDA_FALSE) for x in space}
    expression = substitute_variables_in_expression(expression, substitution)
    return expression.simplify()


def expression_to_space_list(expression: Expression) -> list[dict[str, int]]:
    """
    Convert a PyEDA expression to a list of subspaces whose union represents
    an equivalent set of network states.

    Note that the spaces are not necessarily pair-wise disjoint. Also,
    the list is not necessarily minimal.
    """

    # TODO:
    #  Function `to_dnf` in PyEDA actually produces a minimal
    #  (or at least in some sense canonical) DNF. In the future, we might
    #  want to either enforce this explicitly or relax this requirement.

    sub_spaces: list[dict[str, int]] = []
    expression_dnf = expression.to_dnf()

    for clause in expression_dnf.xs:  # type: ignore
        sub_space: dict[str, int] = {}

        # Since we know this is a DNF clause, it can only be
        # a literal, or a conjunction of literals.
        # TODO: investigate the types here more closely... something strange is going on
        literals = [clause] if isinstance(clause, Literal) else clause.xs  # type: ignore # noqa
        for literal in literals:  # type: ignore
            var = str(literal.inputs[0])  # type: ignore
            if isinstance(literal, Variable):
                sub_space[var] = 1
            elif isinstance(literal, Complement):
                sub_space[var] = 0
            else:
                raise Exception(
                    f"Unreachable: Invalid literal type `{type(literal)}`."  # type: ignore
                )  # type: ignore

        sub_spaces.append(sub_space)

    return sub_spaces


def space_unique_key(space: dict[str, int], network: BooleanNetwork) -> int:
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
        assert var
        var_index: int = var.as_index()
        # Each variable is encoded as two bits, so the total length
        # of the key is 2 * n and the offset of each variable is 2 * index.
        # 00 - unknown; 10 - zero; 11 - one
        key |= (v + 2) << (2 * var_index)
    return key
