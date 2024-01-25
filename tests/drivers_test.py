from balm.drivers import find_single_drivers, find_single_node_LDOIs
from biodivine_aeon import BooleanNetwork, AsynchronousGraph


def test_find_single_node_LDOIs():
    bn = BooleanNetwork.from_bnet(
        """
    S, S
    A, S | B
    B, A
    C, A | D
    D, C
    E, false
    """
    )
    graph = AsynchronousGraph(bn)
    LDOIs = find_single_node_LDOIs(graph)
    assert LDOIs[("S", 0)] == {"S": 0}
    assert LDOIs[("S", 1)] == {"S": 1, "A": 1, "B": 1, "C": 1, "D": 1}
    assert LDOIs[("A", 0)] == {"B": 0}
    assert LDOIs[("A", 1)] == {"A": 1, "B": 1, "C": 1, "D": 1}
    assert LDOIs[("B", 0)] == {}
    assert LDOIs[("B", 1)] == {"B": 1, "A": 1, "C": 1, "D": 1}
    assert LDOIs[("C", 0)] == {"D": 0}
    assert LDOIs[("C", 1)] == {"C": 1, "D": 1}
    assert LDOIs[("D", 0)] == {}
    assert LDOIs[("D", 1)] == {"D": 1, "C": 1}


def test_find_single_drivers():
    bn = BooleanNetwork.from_bnet(
        """
    S, S
    A, S | B
    B, A
    C, A | D
    D, C
    E, true
    """
    )
    graph = AsynchronousGraph(bn)
    LDOIs = find_single_node_LDOIs(graph)
    assert find_single_drivers({"A": 0, "B": 0}, graph) == {("A", 0)}
    assert find_single_drivers({"A": 0, "B": 0}, graph, LDOIs) == {("A", 0)}
    assert find_single_drivers({"A": 1, "B": 1}, graph) == {("A", 1), ("B", 1), ("S", 1)}
    assert find_single_drivers({"A": 1, "B": 1}, graph, LDOIs) == {
        ("A", 1),
        ("B", 1),
        ("S", 1),
    }
    assert find_single_drivers({"C": 0, "D": 0}, graph) == {("C", 0)}
    assert find_single_drivers({"C": 0, "D": 0}, graph, LDOIs) == {("C", 0)}
    assert find_single_drivers({"C": 1, "D": 1}, graph) == {
        ("A", 1),
        ("B", 1),
        ("C", 1),
        ("D", 1),
        ("S", 1),
    }
    assert find_single_drivers({"C": 1, "D": 1}, graph, LDOIs) == {
        ("A", 1),
        ("B", 1),
        ("C", 1),
        ("D", 1),
        ("S", 1),
    }
