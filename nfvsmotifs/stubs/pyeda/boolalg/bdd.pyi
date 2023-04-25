"""
This type stub file was generated by pyright.
"""

from pyeda.boolalg import boolfunc
from pyeda.util import cached_property

"""
The :mod:`pyeda.boolalg.bdd` module implements
Boolean functions represented as binary decision diagrams.

Interface Functions:

* :func:`bddvar` --- Return a unique BDD variable
* :func:`expr2bdd` --- Convert an expression into a binary decision diagram
* :func:`bdd2expr` --- Convert a binary decision diagram into an expression
* :func:`upoint2bddpoint` --- Convert an untyped point into a BDD point
* :func:`ite` --- BDD if-then-else operator

Interface Classes:

* :class:`BDDNode`
* :class:`BinaryDecisionDiagram`

  * :class:`BDDConstant`
  * :class:`BDDVariable`
"""
_VARS = ...
_NODES = ...
_BDDS = ...
class BDDNode:
    """Binary decision diagram node

    A BDD node represents one cofactor in the decomposition of a Boolean
    function.
    Nodes are uniquely identified by a ``root`` integer,
    ``lo`` child node, and ``hi`` child node:

    * ``root`` is the cofactor variable's ``uniqid`` attribute
    * ``lo`` is the zero cofactor node
    * ``hi`` is the one cofactor node

    The ``root`` of the zero node is -1,
    and the ``root`` of the one node is -2.
    Both zero/one nodes have ``lo=None`` and ``hi=None``.

    Do **NOT** create BDD nodes using the ``BDDNode`` constructor.
    BDD node instances are managed internally.
    """
    def __init__(self, root, lo, hi) -> None:
        ...
    


def bddvar(name, index=...): # -> BDDVariable:
    r"""Return a unique BDD variable.

    A Boolean *variable* is an abstract numerical quantity that may assume any
    value in the set :math:`B = \{0, 1\}`.
    The ``bddvar`` function returns a unique Boolean variable instance
    represented by a binary decision diagram.
    Variable instances may be used to symbolically construct larger BDDs.

    A variable is defined by one or more *names*,
    and zero or more *indices*.
    Multiple names establish hierarchical namespaces,
    and multiple indices group several related variables.
    If the ``name`` parameter is a single ``str``,
    it will be converted to ``(name, )``.
    The ``index`` parameter is optional;
    when empty, it will be converted to an empty tuple ``()``.
    If the ``index`` parameter is a single ``int``,
    it will be converted to ``(index, )``.

    Given identical names and indices, the ``bddvar`` function will always
    return the same variable:

    >>> bddvar('a', 0) is bddvar('a', 0)
    True

    To create several single-letter variables:

    >>> a, b, c, d = map(bddvar, 'abcd')

    To create variables with multiple names (inner-most first):

    >>> fifo_push = bddvar(('push', 'fifo'))
    >>> fifo_pop = bddvar(('pop', 'fifo'))

    .. seealso::
       For creating arrays of variables with incremental indices,
       use the :func:`pyeda.boolalg.bfarray.bddvars` function.
    """
    ...

def expr2bdd(expr): # -> BinaryDecisionDiagram:
    """Convert an expression into a binary decision diagram."""
    ...

def bdd2expr(bdd, conj=...): # -> Expression:
    """Convert a binary decision diagram into an expression.

    This function will always return an expression in two-level form.
    If *conj* is ``False``, return a sum of products (SOP).
    Otherwise, return a product of sums (POS).

    For example::

       >>> a, b = map(bddvar, 'ab')
       >>> bdd2expr(~a | b)
       Or(~a, And(a, b))
    """
    ...

def upoint2bddpoint(upoint): # -> dict[Unknown, Unknown]:
    """Convert an untyped point into a BDD point.

    .. seealso::
       For definitions of points and untyped points,
       see the :mod:`pyeda.boolalg.boolfunc` module.
    """
    ...

def ite(f, g, h): # -> BinaryDecisionDiagram:
    r"""BDD if-then-else (ITE) operator

    The *f*, *g*, and *h* arguments are BDDs.

    The ITE(f, g, h) operator means
    "if *f* is true, return *g*, else return *h*".

    It is equivalent to:

    * DNF form: ``f & g | ~f & h``
    * CNF form: ``(~f | g) & (f | h)``
    """
    ...

class BinaryDecisionDiagram(boolfunc.Function):
    """Boolean function represented by a binary decision diagram

    .. seealso::
       This is a subclass of :class:`pyeda.boolalg.boolfunc.Function`

    BDDs have a single attribute, ``node``,
    that points to a node in the managed unique table.

    There are two ways to construct a BDD:

    * Convert an expression using the ``expr2bdd`` function.
    * Use operators on existing BDDs.

    Use the ``bddvar`` function to create BDD variables,
    and use the Python ``~|&^`` operators for NOT, OR, AND, XOR.

    For example::

       >>> a, b, c, d = map(bddvar, 'abcd')
       >>> f = ~a | b & c ^ d

    The ``BinaryDecisionDiagram`` class is useful for type checking,
    e.g. ``isinstance(f, BinaryDecisionDiagram)``.

    Do **NOT** create a BDD using the ``BinaryDecisionDiagram`` constructor.
    BDD instances are managed internally,
    and you will not be able to use the Python ``is`` operator to establish
    formal equivalence with manually constructed BDDs.
    """
    def __init__(self, node) -> None:
        ...
    
    def __invert__(self): # -> BinaryDecisionDiagram:
        ...
    
    def __or__(self, other): # -> BinaryDecisionDiagram:
        ...
    
    def __and__(self, other): # -> BinaryDecisionDiagram:
        ...
    
    def __xor__(self, other): # -> BinaryDecisionDiagram:
        ...
    
    def __rshift__(self, other): # -> BinaryDecisionDiagram:
        ...
    
    def __rrshift__(self, other): # -> BinaryDecisionDiagram:
        ...
    
    @cached_property
    def support(self): # -> frozenset[Any]:
        ...
    
    @cached_property
    def inputs(self): # -> tuple[Unknown, ...]:
        ...
    
    def restrict(self, point): # -> BinaryDecisionDiagram:
        ...
    
    def compose(self, mapping): # -> BinaryDecisionDiagram:
        ...
    
    def satisfy_one(self): # -> dict[Unknown, int] | None:
        ...
    
    def satisfy_all(self): # -> Generator[dict[Unknown, int], None, None]:
        ...
    
    def is_zero(self): # -> bool:
        ...
    
    def is_one(self): # -> bool:
        ...
    
    @staticmethod
    def box(obj): # -> BinaryDecisionDiagram | BDDConstant:
        ...
    
    def dfs_preorder(self): # -> Generator[Unknown, None, None]:
        """Iterate through nodes in depth first search (DFS) pre-order."""
        ...
    
    def dfs_postorder(self): # -> Generator[Unknown, None, None]:
        """Iterate through nodes in depth first search (DFS) post-order."""
        ...
    
    def bfs(self): # -> Generator[Unknown, None, None]:
        """Iterate through nodes in breadth first search (BFS) order."""
        ...
    
    def equivalent(self, other): # -> bool:
        """Return whether this BDD is equivalent to *other*.

        You can also use Python's ``is`` operator for BDD equivalency testing.

        For example::

           >>> a, b, c = map(bddvar, 'abc')
           >>> f1 = a ^ b ^ c
           >>> f2 = a & ~b & ~c | ~a & b & ~c | ~a & ~b & c | a & b & c
           >>> f1 is f2
           True
        """
        ...
    
    def to_dot(self, name=...): # -> str:
        """Convert to DOT language representation.

        See the
        `DOT language reference <http://www.graphviz.org/content/dot-language>`_
        for details.
        """
        ...
    


class BDDConstant(BinaryDecisionDiagram):
    """Binary decision diagram constant zero/one

    The ``BDDConstant`` class is useful for type checking,
    e.g. ``isinstance(f, BDDConstant)``.

    Do **NOT** create a BDD using the ``BDDConstant`` constructor.
    BDD instances are managed internally,
    and the BDD zero/one instances are singletons.
    """
    def __init__(self, node, value) -> None:
        ...
    
    def __bool__(self): # -> bool:
        ...
    
    def __int__(self) -> int:
        ...
    
    def __repr__(self): # -> str:
        ...
    
    def __str__(self) -> str:
        ...
    


class BDDVariable(boolfunc.Variable, BinaryDecisionDiagram):
    """Binary decision diagram variable

    The ``BDDVariable`` class is useful for type checking,
    e.g. ``isinstance(f, BDDVariable)``.

    Do **NOT** create a BDD using the ``BDDVariable`` constructor.
    Use the :func:`bddvar` function instead.
    """
    def __init__(self, bvar) -> None:
        ...
    


