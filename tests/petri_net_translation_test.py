from nfvsmotifs.petri_net_translation import sanitize_network_names, network_to_petrinet
from biodivine_aeon import RegulatoryGraph, BooleanNetwork # type: ignore
from networkx import DiGraph, is_isomorphic # type: ignore
import pytest

def test_sanitization():
    rg = RegulatoryGraph([r"a_45[x]", r"b12{z}", "c[", "c]"])
    bn = BooleanNetwork(rg)
    a = bn.variables()[0]
    b = bn.variables()[1]
    c1 = bn.variables()[2]
    c2 = bn.variables()[3]
    assert bn.get_variable_name(a) == r"a_45[x]"
    assert bn.get_variable_name(b) == r"b12{z}"
    assert bn.get_variable_name(c1) == "c["
    assert bn.get_variable_name(c2) == "c]"

    bn = sanitize_network_names(bn)
    assert bn.get_variable_name(a) == "a_45_x_"
    assert bn.get_variable_name(b) == "b12_z_"
    # Name clash is resolved by adding an ID.
    assert bn.get_variable_name(c1) == "c_"
    assert bn.get_variable_name(c2) == "c__id4"

def test_sanitization_failing():
    rg = RegulatoryGraph(["x_", "x_id2", "x["])
    bn = BooleanNetwork(rg)
    try:
        bn = sanitize_network_names(bn)
        pytest.fail("This network should not be sanitizeable.")
    except:
        pass

def test_translation():
    # A very very simple network for which we know how the 
    # translation should look like.
    bn = BooleanNetwork.from_bnet("""
        A, !A & B
        B, !B & !A
    """)

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

    assert is_isomorphic(network_to_petrinet(bn), expected)