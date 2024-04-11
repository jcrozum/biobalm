# type: ignore
import pytest
from biodivine_aeon import BooleanNetwork
from networkx import DiGraph, is_isomorphic  # type: ignore

from biobalm.petri_net_translation import (
    extract_variable_names,
    network_to_petrinet,
    sanitize_network_names,
)


def test_sanitization():
    bn = BooleanNetwork([r"a_45[x]", r"b12{z}", "c[", "c]"])
    a, b, c1, c2 = bn.variables()
    assert bn.variable_names() == [r"a_45[x]", r"b12{z}", "c[", "c]"]

    bn = sanitize_network_names(bn)
    assert bn.variable_names() == ["a_45_x_", "b12_z_", "c_", "_c_"]


def test_sanitization_failing():
    bn = BooleanNetwork(["x_", "x_id2", "x["])
    try:
        bn = sanitize_network_names(bn)
        pytest.fail("This network should not be sanitizeable.")
    except:  # noqa: E722
        pass


def test_translation():
    # A very very simple network for which we know how the
    # translation should look like.
    bn = BooleanNetwork.from_bnet(
        """
        A, !A & B
        B, !B & !A
    """
    )

    expected = DiGraph()
    expected.add_node("b0_A", kind="place")
    expected.add_node("b1_A", kind="place")
    expected.add_node("b0_B", kind="place")
    expected.add_node("b1_B", kind="place")

    # A goes up if A=0 and B=1.
    expected.add_node("tr_A_up_1", kind="transition")
    expected.add_edge("b0_A", "tr_A_up_1")
    expected.add_edge("tr_A_up_1", "b1_A")
    expected.add_edge("b1_B", "tr_A_up_1")
    expected.add_edge("tr_A_up_1", "b1_B")

    # A goes down if A=1, regardless of B.
    expected.add_node("tr_A_down_1", kind="transition")
    expected.add_edge("b1_A", "tr_A_down_1")
    expected.add_edge("tr_A_down_1", "b0_A")

    # B goes up if B=0 and A=0.
    expected.add_node("tr_B_up_1", kind="transition")
    expected.add_edge("b0_B", "tr_B_up_1")
    expected.add_edge("tr_B_up_1", "b1_B")
    expected.add_edge("b0_A", "tr_B_up_1")
    expected.add_edge("tr_B_up_1", "b0_A")

    # B goes down if B=1, regardless of A.
    expected.add_node("tr_B_down_1", kind="transition")
    expected.add_edge("b1_B", "tr_B_down_1")
    expected.add_edge("tr_B_down_1", "b0_B")

    pn = network_to_petrinet(bn)
    assert ["A", "B"] == extract_variable_names(pn)
    assert is_isomorphic(pn, expected)
