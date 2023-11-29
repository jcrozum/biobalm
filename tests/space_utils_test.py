from typing import cast

from biodivine_aeon import BooleanNetwork  # type: ignore
from pyeda.boolalg.expr import Expression, expr

from nfvsmotifs.space_utils import (
    expression_to_space_list,
    is_subspace,
    is_syntactic_trap_space,
    percolate_network,
    percolate_pyeda_expression,
    percolate_space,
    space_unique_key,
)


def test_is_subspace():
    assert is_subspace({"x": 0, "y": 1}, {"x": 0})
    assert not is_subspace({"x": 1, "y": 0}, {"x": 0, "y": 0})
    assert is_subspace({"x": 0, "y": 1}, {"x": 0})
    assert not is_subspace({"x": 1, "y": 0}, {"x": 0, "y": 0})


def test_expression_percolation():
    e = cast(Expression, expr("(a & ~x) | (a & y)"))

    assert expr("a") == percolate_pyeda_expression(e, {"x": 0, "y": 1})
    assert expr(False) == percolate_pyeda_expression(e, {"a": 0})
    assert expr("a") == percolate_pyeda_expression(e, {"x": 0, "y": 1})
    assert expr(False) == percolate_pyeda_expression(e, {"a": 0})


def test_space_percolation():
    bn = BooleanNetwork.from_bnet(
        """
        a, b
        b, c
        c, a
    """
    )

    assert {"a": 0, "b": 0, "c": 0} == percolate_space(bn, {"a": 0})
    # assert {} == percolate_space(bn, {"a": 0})[1]
    assert is_syntactic_trap_space(bn, {"a": 0, "b": 0, "c": 0})
    assert {"a": 1, "b": 1, "c": 1} == percolate_space(bn, {"a": 1})
    # assert {} == percolate_space(bn, {"a": 1})[1]
    assert is_syntactic_trap_space(bn, {"a": 1, "b": 1, "c": 1})

    bn = BooleanNetwork.from_bnet(
        """
    a, b
    b, !c
    c, a
    """
    )

    assert {"a": 0, "b": 0, "c": 0} == percolate_space(
        bn, {"a": 0, "b": 0, "c": 0}, strict_percolation=False
    )
    assert {"a": 0, "c": 0} == percolate_space(
        bn, {"a": 0, "b": 0, "c": 0}, strict_percolation=True
    )

    # The conflict is on b. The rest is fine.
    # assert {"b": 1} == percolate_space(bn, {"a": 0, "b": 0, "c": 0})[1]
    assert not is_syntactic_trap_space(bn, {"a": 0})
    assert is_syntactic_trap_space(bn, {})

    bn = BooleanNetwork.from_bnet(
        """
    a, !b
    b, a
    c, a & c & d | b & !c | c & !d
    d, !a | d
    """
    )
    assert {"b": 1, "c": 1} == percolate_space(bn, {"a": 1})


def test_constant_percolation():
    bn = BooleanNetwork.from_bnet(
        """
        a, true
        b, b
        c, a | b
    """
    )

    assert {"a": 1, "c": 1} == percolate_space(bn, {}, strict_percolation=False)
    # assert {"a": 1} == percolate_space(bn, {"a": 0}, strict_percolation=False)[1]
    assert {} == percolate_space(bn, {}, strict_percolation=True)
    # assert {} == percolate_space(bn, {"a": 0}, strict_percolation=True)[1]


def test_network_percolation():
    bn = BooleanNetwork.from_bnet(
        """
        a, c & b
        b, !a
        c, c
    """
    )

    percolated_bn = percolate_network(bn, percolate_space(bn, {"c": 0}))
    percolated_bn = percolate_network(bn, percolate_space(bn, {"c": 0}))
    assert "false" == percolated_bn.get_update_function("c")
    assert "false" == percolated_bn.get_update_function("a")
    assert "true" == percolated_bn.get_update_function("b")
    percolated_bn = percolate_network(bn, percolate_space(bn, {"c": 1}))
    percolated_bn = percolate_network(bn, percolate_space(bn, {"c": 1}))
    assert "true" == percolated_bn.get_update_function("c")
    assert "b" == percolated_bn.get_update_function("a")
    assert "!a" == percolated_bn.get_update_function("b")


def test_expression_to_spaces():
    e = cast(Expression, expr("(a & c) | (~d & (a | c)) | f"))

    spaces = expression_to_space_list(e)

    assert {"f": 1} in spaces
    assert {"a": 1, "c": 1} in spaces
    assert {"d": 0, "a": 1} in spaces
    assert {"d": 0, "c": 1} in spaces
    assert {"a": 1, "c": 0} not in spaces


def test_space_unique_key():
    bn = BooleanNetwork.from_bnet(
        """
        a, c & b
        b, !a
        c, c
    """
    )
    assert space_unique_key({"a": 1}, bn) == space_unique_key({"a": 1}, bn)
    assert space_unique_key({"a": 1}, bn) != space_unique_key({"b": 1}, bn)
