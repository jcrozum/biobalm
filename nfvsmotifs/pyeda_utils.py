from __future__ import annotations
"""
    Some utility methods, mainly for converting and modifying PyEDA expressions.
"""

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from pyeda.inter import Expression # type: ignore

import pyeda.boolalg.expr as pyeda_expression # type: ignore
from pyeda.inter import Not, And, Or, Equal, Xor, Implies # type: ignore

PYEDA_TRUE = pyeda_expression.expr(1)
PYEDA_FALSE = pyeda_expression.expr(0)

def substitute_variables_in_expression(expression: Expression, items: dict[str, Expression]) -> Expression:    
    """
        A substitution method which replaces all occurences of the variables specified in `items`
        with the corresponding `Expression` objects.
    """
    if expression == PYEDA_TRUE or expression == PYEDA_FALSE:
        # Keep constants.
        return expression
    if type(expression) == pyeda_expression.Variable:
        # Positive literals are resolved through `items` if possible.
        key = str(expression)
        if key in items:
            return items[key]
        else:
            return expression
    if type(expression) == pyeda_expression.Complement:
        # Complement is just a negative literal.
        key = str(expression.inputs[0])
        if key in items:
            return Not(items[key])
        else:
            return expression
    if type(expression) == pyeda_expression.NotOp:
        inner = substitute_variables_in_expression(expression.x, items)
        return Not(expression)
    if type(expression) == pyeda_expression.AndOp:
        inner = [substitute_variables_in_expression(x, items) for x in expression.xs]
        return And(*inner)
    if type(expression) == pyeda_expression.OrOp:
        inner = [substitute_variables_in_expression(x, items) for x in expression.xs]
        return Or(*inner)
    if type(expression) == pyeda_expression.EqualOp:
        inner = [substitute_variables_in_expression(x, items) for x in expression.xs]
        return Equal(*inner)
    if type(expression) == pyeda_expression.XorOp:
        inner = [substitute_variables_in_expression(x, items) for x in expression.xs]
        return Xor(inner)
    if type(expression) == pyeda_expression.ImpliesOp:
        p = substitute_variables_in_expression(expression.xs[0], items)
        q = substitute_variables_in_expression(expression.xs[1], items)
        return Implies(p, q)
    raise Exception(f"Unknown PyEDA operator: {type(expression)}.")

def pyeda_to_aeon(expression: Expression) -> str:
    """
        Convert a PyEDA expression to a string representation that is valid in AEON.py 
        (but also other formats, like .bnet).

        Note: Right now, I have not found a better way to do this. If you find something
        simpler, go for it...
    """
    if expression == PYEDA_TRUE:
        return "true"
    if expression == PYEDA_FALSE:        
        return "false"
    if type(expression) == pyeda_expression.Variable:
        return str(expression)        
    if type(expression) == pyeda_expression.Complement:
        return f"!{str(expression.inputs[0])}"        
    if type(expression) == pyeda_expression.NotOp:
        return f"!{pyeda_to_aeon(expression.x)}"        
    if type(expression) == pyeda_expression.AndOp:
        if len(expression.xs) == 0:
            return "true"
        inner = " & ".join([pyeda_to_aeon(x) for x in expression.xs])
        return f"({inner})"
    if type(expression) == pyeda_expression.OrOp:
        if len(expression.xs) == 0:
            return "false"
        inner = " | ".join([pyeda_to_aeon(x) for x in expression.xs])
        return f"({inner})"
    if type(expression) == pyeda_expression.EqualOp:
        if len(expression.xs) == 0:
            return "true"
        inner = " <=> ".join([pyeda_to_aeon(x) for x in expression.xs])
        return f"({inner})"
    if type(expression) == pyeda_expression.XorOp:
        if len(expression.xs) == 0:
            return "false"
        inner = " ^ ".join([pyeda_to_aeon(x) for x in expression.xs])
        return f"({inner})"
    if type(expression) == pyeda_expression.ImpliesOp:
        p = pyeda_to_aeon(expression.xs[0])
        q = pyeda_to_aeon(expression.xs[1])
        return f"({p} => {q})"
    raise Exception(f"Unknown PyEDA operator: {type(expression)}.")

def aeon_to_pyeda(expression: str) -> Expression:
    """
        Convert a Boolean expression from AEON.py to PyEDA.
    """    
    return pyeda_expression.expr(expression.replace("!", "~"))