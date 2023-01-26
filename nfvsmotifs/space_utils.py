from __future__ import annotations
"""
    Some basic utility operations on spaces (partial assignments of BN variables).

    Each space is represented as a dictionary with a subset of variable names as 
    keys and values `'0'`/`'1'` assigned to fixed variables.
"""

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from pyeda.inter import Expression # type: ignore

from pyeda.inter import expr # type: ignore
from biodivine_aeon import BooleanNetwork, RegulatoryGraph # type: ignore
from nfvsmotifs.pyeda_utils import aeon_to_pyeda

from nfvsmotifs.pyeda_utils import substitute_variables_in_expression, pyeda_to_aeon, aeon_to_pyeda, PYEDA_TRUE, PYEDA_FALSE

def is_subspace(x: dict[str, str], y: dict[str, str]) -> bool:
    """
        Checks if `x` is a subspace of `y`.
    """
    for var in y:
        if var not in x:
            return False
        if x[var] != y[var]:
            return False
    return True

def is_syntactic_trap_space(bn: BooleanNetwork, space: dict[str, str]) -> bool:
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
            if space[var_name] != str(expression):
                print(space[var_name], str(expression), bn.get_update_function(var), space)
                return False
    return True

def percolate_space(network: BooleanNetwork, space: dict[str, str]) -> tuple[dict[str, str], dict[str, str]]:
    """
        Takes a Boolean network and a space (partial assignment of `"0"`/`"1"` 
        to the network variables). It then percolates the values in the given 
        `space` to the remaining network variables based on the update functions
        of the given `network`. 
        
        If the argument is a trap space, then the result is a subspace of 
        the argument and is also a trap space. 
        
        However, when the argument is a general space, the percolation can 
        actually lead "outside" of the original space. In such case, the original 
        fixed value is *not* modified and the conflict will remain in the 
        resulting space.

        We then return these percolated values for the conflicting variables
        as a second member of the result tuple.
    """
    # Make a copy of the original space
    result = { var:space[var] for var in space }
    conflicts = {}
    done = False
    while not done:
        done = True
        for var in network.variables():
            var_name = network.get_variable_name(var)            
            expression = aeon_to_pyeda(network.get_update_function(var))
            
            # if the var is already constant, it doesn't count
            if expression == PYEDA_TRUE or expression == PYEDA_FALSE: 
                continue
            
            expression = percolate_pyeda_expression(expression, result)
            if expression == PYEDA_TRUE or expression == PYEDA_FALSE:                
                if var_name not in result:
                    # Fortunately, PyEDA resolves true as '1' and false as '0', 
                    # so we can use this directly.
                    result[var_name] = str(expression)
                    done = False
                if var_name in result and result[var_name] != str(expression):
                    conflicts[var_name] = str(expression)
    
    return (result, conflicts)

def percolate_network(bn: BooleanNetwork, space: dict[str, str]) -> BooleanNetwork:
    """
        Takes an AEON.py Boolean network and a space (partial assignment of
        network variables to `'0'`/`'1'`). It then produces a new network with
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
        reg['observable'] = False
        if 'monotonicity' in reg:
            del reg['monotonicity']
        new_rg.add_regulation(reg)
    
    # Copy the Boolean network, but with simplified expressions.
    new_bn = BooleanNetwork(new_rg)
    for var in bn.variables():
        name = bn.get_variable_name(var)
        new_var = new_bn.find_variable(name)    # The ids should be the same, but just in case.
        new_expr = None
        if name in space:
            # If the value is fixed, just use it as a value directly.
            if space[name] == "1":
                new_expr = PYEDA_TRUE
            if space[name] == "0":
                new_expr = PYEDA_FALSE
        else:
            # If the value is not fixed, use a simplified expression.
            expression = aeon_to_pyeda(bn.get_update_function(var))
            new_expr = percolate_pyeda_expression(expression, space)

        new_bn.set_update_function(var, pyeda_to_aeon(new_expr))        

    return new_bn

def percolate_pyeda_expression(expression: Expression, space: dict[str, str]) -> Expression:
    """
        Takes a PyEDA expression and a subspace (dictionary assigning `"1"`/`"0"` to
        a subset of variables). Returns a simplified expression that is valid
        for exactly the same members of the given `space` as the original expression. 
        The resulting expression does not depend on the variables which are fixed 
        in the given `space`.
    """
    substitution = { x: expr(space[x]) for x in space }
    expression = substitute_variables_in_expression(expression, substitution)
    return expression.simplify()
