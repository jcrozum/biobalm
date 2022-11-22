
from biodivine_aeon import BooleanNetwork # type:ignore
from interaction_graph_utils import infer_signed_interaction_graph

def test_ig_inference():
    bn = BooleanNetwork.from_bnet("""
        # Just a normal function.
        b, a | !b
        # Contradiciton on `a` - the regulation should not appear in the result
        # Also, non-monotonic dependence on b and c.
        a, (a & !a) | (b <=> c)
        c, c
    """)
    ig = infer_signed_interaction_graph(bn)

    edges = { edge:ig.get_edge_data(edge[0], edge[1])['sign'] for edge in ig.edges }
    assert len(edges) == 5
    assert edges[('a', 'b')] == "+"
    assert edges[('b', 'b')] == "-"
    assert edges[('b', 'a')] == "?"
    assert edges[('c', 'a')] == "?"
    assert edges[('c', 'c')] == "+"
    assert ('a', 'a') not in edges
