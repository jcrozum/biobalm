from biodivine_aeon import BooleanNetwork, BooleanExpression, AsynchronousGraph

from balm.space_utils import (
    expression_to_space_list,
    is_subspace,
    percolate_network,
    percolation_conflicts,
    percolate_expression,
    percolate_space,
    percolate_space_strict,
    space_unique_key,
)


def test_is_subspace():
    assert is_subspace({"x": 0, "y": 1}, {"x": 0})
    assert not is_subspace({"x": 1, "y": 0}, {"x": 0, "y": 0})
    assert is_subspace({"x": 0, "y": 1}, {"x": 0})
    assert not is_subspace({"x": 1, "y": 0}, {"x": 0, "y": 0})


def test_expression_percolation():
    e = BooleanExpression("(a & !x) | (a & y)")

    assert BooleanExpression("a") == percolate_expression(e, {"x": 0, "y": 1})
    assert BooleanExpression("false") == percolate_expression(e, {"a": 0})
    assert BooleanExpression("a") == percolate_expression(e, {"x": 0, "y": 1})
    assert BooleanExpression("false") == percolate_expression(e, {"a": 0})


def test_space_percolation():
    bn = BooleanNetwork.from_bnet(
        """
        a, b
        b, c
        c, a
    """
    )
    graph = AsynchronousGraph(bn)

    assert {"a": 0, "b": 0, "c": 0} == percolate_space(graph, {"a": 0})
    assert {} == percolation_conflicts(graph, {"a": 0})
    assert {"a": 1, "b": 1, "c": 1} == percolate_space(graph, {"a": 1})
    assert {} == percolation_conflicts(graph, {"a": 1})
    
    bn = BooleanNetwork.from_bnet(
        """
    a, b
    b, !c
    c, a
    """
    )
    graph = AsynchronousGraph(bn)

    assert {"a": 0, "b": 0, "c": 0} == percolate_space(graph, {"a": 0, "b": 0, "c": 0})
    assert {"a": 0, "c": 0} == percolate_space_strict(graph, {"a": 0, "b": 0, "c": 0})

    # The conflict is on b. The rest is fine.
    assert {"b": 1} == percolation_conflicts(graph, {"a": 0, "b": 0, "c": 0})
    
    bn = BooleanNetwork.from_bnet(
        """
    a, !b
    b, a
    c, a & c & d | b & !c | c & !d
    d, !a | d
    """
    )
    graph = AsynchronousGraph(bn)
    assert {"b": 1, "c": 1} == percolate_space_strict(graph, {"a": 1})


def test_constant_percolation():
    bn = BooleanNetwork.from_bnet(
        """
        a, true
        b, b
        c, a | b
    """
    )
    graph = AsynchronousGraph(bn)

    assert {"a": 1, "c": 1} == percolate_space(graph, {})
    assert {"a": 1} == percolation_conflicts(graph, {"a": 0}, strict_percolation=False)
    assert {} == percolate_space_strict(graph, {})
    assert {} == percolation_conflicts(graph, {"a": 0}, strict_percolation=True)


def test_network_percolation():
    bn = BooleanNetwork.from_bnet(
        """
        a, c & b
        b, !a
        c, c
    """
    )
    graph = AsynchronousGraph(bn)

    percolated_bn = percolate_network(bn, percolate_space(graph, {"c": 0}), ctx=graph.symbolic_context())
    percolated_bn = percolate_network(bn, percolate_space(graph, {"c": 0}), ctx=graph.symbolic_context())
    assert "false" == str(percolated_bn.get_update_function("c"))
    assert "false" == str(percolated_bn.get_update_function("a"))
    assert "true" == str(percolated_bn.get_update_function("b"))
    percolated_bn = percolate_network(bn, percolate_space(graph, {"c": 1}), ctx=graph.symbolic_context())
    percolated_bn = percolate_network(bn, percolate_space(graph, {"c": 1}), ctx=graph.symbolic_context())
    assert "true" == str(percolated_bn.get_update_function("c"))
    assert "b" == str(percolated_bn.get_update_function("a"))
    assert "!a" == str(percolated_bn.get_update_function("b"))


def test_expression_to_spaces():
    e = BooleanExpression("(a & c) | (!d & (a | c)) | f")

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
