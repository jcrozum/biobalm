from biodivine_aeon import AsynchronousGraph, BooleanExpression, BooleanNetwork

from biobalm.space_utils import (
    expression_to_space_list,
    is_subspace,
    percolate_network,
    percolate_space,
    percolate_space_strict,
    percolation_conflicts,
    restrict_expression,
    space_unique_key,
)


def test_is_subspace():
    assert is_subspace({"x": 0, "y": 1}, {"x": 0})
    assert not is_subspace({"x": 1, "y": 0}, {"x": 0, "y": 0})
    assert is_subspace({"x": 0, "y": 1}, {"x": 0})
    assert not is_subspace({"x": 1, "y": 0}, {"x": 0, "y": 0})


def test_expression_percolation():
    e = BooleanExpression("(a & !x) | (a & y)")

    assert "a" == str(restrict_expression(e, {"x": 0, "y": 1}))
    assert "false" == str(restrict_expression(e, {"a": 0}))
    assert "a" == str(restrict_expression(e, {"x": 0, "y": 1}))
    assert "false" == str(restrict_expression(e, {"a": 0}))


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
    assert set() == percolation_conflicts(graph, {"a": 0})
    assert {"a": 1, "b": 1, "c": 1} == percolate_space(graph, {"a": 1})
    assert set() == percolation_conflicts(graph, {"a": 1})

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
    assert {"b"} == percolation_conflicts(
        graph, {"a": 0, "b": 0, "c": 0}, strict_percolation=False
    )

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
    assert {"a"} == percolation_conflicts(graph, {"a": 0}, strict_percolation=False)
    assert {} == percolate_space_strict(graph, {})
    assert set() == percolation_conflicts(graph, {"a": 0}, strict_percolation=True)


def test_network_percolation():
    bn = BooleanNetwork.from_bnet(
        """
        a, c & b
        b, !a
        c, c
    """
    )
    graph = AsynchronousGraph(bn)

    percolated_bn = percolate_network(
        bn, percolate_space(graph, {"c": 0}), symbolic_network=graph
    )
    percolated_bn = percolate_network(
        bn, percolate_space(graph, {"c": 0}), symbolic_network=graph
    )
    assert "false" == str(percolated_bn.get_update_function("c"))
    assert "false" == str(percolated_bn.get_update_function("a"))
    assert "true" == str(percolated_bn.get_update_function("b"))
    percolated_bn = percolate_network(
        bn, percolate_space(graph, {"c": 1}), symbolic_network=graph
    )
    percolated_bn = percolate_network(
        bn, percolate_space(graph, {"c": 1}), symbolic_network=graph
    )
    assert "true" == str(percolated_bn.get_update_function("c"))
    assert "b" == str(percolated_bn.get_update_function("a"))
    assert "!a" == str(percolated_bn.get_update_function("b"))


def test_expression_to_spaces():
    e = BooleanExpression("(a & c) | (!d & (a | c)) | f")

    spaces = expression_to_space_list(e)

    assert {"a": 0, "c": 0, "f": 1} in spaces
    assert {"a": 0, "c": 1, "d": 0} in spaces
    assert {"a": 0, "c": 1, "d": 1, "f": 1} in spaces
    assert {"a": 1, "c": 0, "d": 0} in spaces
    assert {"a": 1, "c": 0, "d": 1, "f": 1} in spaces
    assert {"a": 1, "c": 1} in spaces

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


def test_perc_and_remove_constants_from_bn():
    bn = BooleanNetwork.from_bnet(
        """targets,factors
    constant1_1, (constant1_1 | !constant1_1)
    constant1_0, (constant1_0 & !constant1_0)
    constant2_1, true
    constant2_0, false
    source, source
    oscillator, !oscillator
    source_after_perc, source_after_perc & constant1_1
    after_perc_0, after_perc_0 & constant1_0"""
    ).infer_valid_graph()

    clean_bnet = percolate_network(bn, {}, remove_constants=True).to_bnet()

    bn2 = BooleanNetwork.from_bnet(
        """targets,factors
    source, source
    oscillator, !oscillator
    source_after_perc, source_after_perc"""
    )

    clean_bnet2 = percolate_network(bn2, {}, remove_constants=True).to_bnet()

    assert clean_bnet == clean_bnet2
