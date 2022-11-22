"""
    Some basic utility operations on spaces (partial assignments of BN variables).

    Each space is represented as a dictionary with a subset of variable names as 
    keys and values `'0'`/`'1'` assigned to fixed variables.
"""

import pyeda.boolalg.expr as pyeda_expression
from pyeda.inter import *
from biodivine_aeon import BooleanNetwork, RegulatoryGraph

def is_subspace(x, y) -> bool:
    """
        Checks if `x` is a subspace of `y`.
    """
    for var in y:
        if var not in x:
            return False
        if x[var] != y[var]:
            return False
    return True

def is_trap_space(bn: BooleanNetwork, space) -> bool:
    """
        Checks if the given `space` is a trap space in the given `BooleanNetwork`.
    """
    for var in bn.variables():
        var_name = bn.get_variable_name(var)

        if var_name in space:
            expression = expr(bn.get_update_function(var).replace("!", "~"))
            expression = percolate_pyeda_expression(expression, space)
            if space[var_name] != str(expression):
                return False
    return True

def percolate_space(network: BooleanNetwork, space: dict) -> dict:
    """
        Takes a Boolean network and a space (partial assignment of `"0"`/`"1"` 
        to the network variables). It then percolates the values in the given 
        `space` to the remaining network variables based on the update functions
        of the given `network`. 
        
        If the argument is a trap space, then the result is a subspace of 
        the argument and is also a trap space. However, when the argument is
        a general space, the result is still a subspace, but not necessarily
        a trap.

        Additionally, if the argument is not a trap space, the result
        can percolate completely outside of the given space, in which case
        we return `None`.
    """
    # Make a copy of the original space
    result = { var:space[var] for var in space }
    done = False
    while not done:
        done = True
        for var in network.variables():
            var_name = network.get_variable_name(var)            
            expression = expr(network.get_update_function(var).replace("!", "~"))
            expression = percolate_pyeda_expression(expression, result)
            is_constant = expression == expr(True) or expression == expr(False)
            if expression == expr(True) or expression == expr(False):                
                if var_name not in result:
                    # Fortunately, PyEDA resolved true as '1' and false as '0', 
                    # so we can use this directly.
                    result[var_name] = str(expression)
                    done = False
                elif result[var_name] != str(expression):
                    # This space percolates completely outside of the original space.
                    return None
    return result

def percolate_network(bn: BooleanNetwork, space: dict) -> BooleanNetwork:
    """
        Takes an AEON.py Boolean network and a space (partial assignment of
        network variables to `'0'`/`'1'`). It then produces a new network with
        update functions percolated based on the supplied space.

        There are two caveats to this operation:
        
            (1) The given space must be a trap space. Otherwise, there would
            be no meaningful relationship between the dynamics of the original
            and the percolated network. If the argument is not a trap space,
            the function returns `None`.

            (2) The underlying regulatory graph of the new network retains all 
            regulations of the original network, but all integrity constraints 
            (essentiality, monotonicity) are removed, because they most likely 
            no longer hold in the new network.
    """

    percolated = percolate_space(bn, space)
    if not (percolated and is_trap_space(bn, percolated)):
        return None

    # Use the percolated space, just in case it is smaller.
    space = percolated
    
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
        expression = expr(bn.get_update_function(var).replace("!", "~"))
        expression = percolate_pyeda_expression(expression, space)
        # True/False and negation in PyEDA are different, rest should be the same.
        if expression == expr(True):
            new_bn.set_update_function(new_var, "true")
        elif expression == expr(False):
            new_bn.set_update_function(new_var, "false")
        else:
            new_bn.set_update_function(new_var, str(expression).replace("~", "!"))        

    return new_bn

def percolate_pyeda_expression(expression, space: dict):
    """
        Takes a PyEDA expression and a subspace (dictionary assigning `"1"`/`"0"` to
        a subset of variables). Returns a simplified expression that is valid
        for exactly the same members of the given `space` as the original expression. 
        The resulting expression does not depend on the variables which are fixed 
        in the given `space`.
    """

    def substitute(expression, space: dict):
        """
            Substitutes the know constants from the given space into the PyEDA expression.
        """
        if type(expression) == pyeda_expression._One or type(expression) == pyeda_expression._Zero:
            # Keep constants.
            return expression
        if type(expression) == pyeda_expression.Variable:
            # Positive literals are resolved through `space` if possible.
            key = str(expression)
            if key in space:
                return expr(space[key])                
            else:
                return expression
        if type(expression) == pyeda_expression.Complement:
            # Complement is just a negative literal.
            key = str(expression.inputs[0])
            if key in space:
                return Not(expr(space[key]))                
            else:
                return expression
        if type(expression) == pyeda_expression.NotOp:
            inner = substitute(expression.x, space)
            return Not(expression)
        if type(expression) == pyeda_expression.AndOp:
            inner = [substitute(x, space) for x in expression.xs]
            return And(*inner)
        if type(expression) == pyeda_expression.OrOp:
            inner = [substitute(x, space) for x in expression.xs]
            return Or(*inner)
        if type(expression) == pyeda_expression.EqualOp:
            inner = [substitute(x, space) for x in expression.xs]
            return Equal(*inner)
        if type(expression) == pyeda_expression.XorOp:
            inner = [substitute(x, space) for x in expression.xs]
            return Xor(inner)
        if type(expression) == pyeda_expression.ImpliesOp:
            p = substitute(expression.xs[0], space)
            q = substitute(expression.xs[1], space)
            return Implies(p, q)
        raise Exception(f"Unknown PyEDA operator: {type(expression)}.")
	
    expression = substitute(expression, space)
    return expression.simplify()


if __name__ == '__main__':
    assert is_subspace({'x': '0', 'y': '1'}, {'x': '0'})
    assert not is_subspace({'x': '1', 'y': '0'}, {'x':'0', 'y':'0'})

    space = {'x': '0', 'y': '1'}
    e = expr("(a & ~x) | (a & y)")
    print(e)

    print(percolate_pyeda_expression(e, {'x': '0', 'y': '1'}))
    assert expr("a") == percolate_pyeda_expression(e, {'x': '0', 'y': '1'})
    print(percolate_pyeda_expression(e, {'a': '0'}))
    assert expr(False) == percolate_pyeda_expression(e, {'a': '0'})

    bn = BooleanNetwork.from_bnet("""
    a, b
    b, c
    c, a
    """)

    print(percolate_space(bn, {'a': '0'}))
    assert {'a': '0', 'b': '0', 'c': '0'} == percolate_space(bn, {'a': '0'})
    assert is_trap_space(bn, {'a': '0', 'b': '0', 'c': '0'})
    print(percolate_space(bn, {'a': '1'}))
    assert {'a': '1', 'b': '1', 'c': '1'} == percolate_space(bn, {'a': '1'})
    assert is_trap_space(bn, {'a': '1', 'b': '1', 'c': '1'})

    bn = BooleanNetwork.from_bnet("""
    a, b
    b, !c
    c, a
    """)

    print(percolate_space(bn, {'a': '0', 'b': '0', 'c': '0'}))
    assert None == percolate_space(bn, {'a': '0', 'b': '0', 'c': '0'})
    assert not is_trap_space(bn, {'a': '0'})
    assert is_trap_space(bn, {})

    bn = BooleanNetwork.from_bnet("""
    a, c & b
    b, !a
    c, c
    """)
    percolated_bn = percolate_network(bn, {'c': '0'})
    assert "false" == percolated_bn.get_update_function("c")
    assert "false" == percolated_bn.get_update_function("a")
    assert "true" == percolated_bn.get_update_function("b")
    percolated_bn = percolate_network(bn, {'c': '1'})
    assert "true" == percolated_bn.get_update_function("c")
    assert "b" == percolated_bn.get_update_function("a")
    assert "!a" == percolated_bn.get_update_function("b")

    print("All checks succeeded.")