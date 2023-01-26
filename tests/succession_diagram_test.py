from nfvsmotifs.SuccessionDiagram import SuccessionDiagram
from biodivine_aeon import BooleanNetwork

def test_succession_diagram_structure():
    bn = BooleanNetwork.from_bnet("""
        x1, x2
        x2, x1
        x3, !x3
        """)

    SD = SuccessionDiagram(bn)
    assert SD.G.number_of_nodes() == 3
    assert SD.G.number_of_edges() == 2
    assert max(d['depth'] for n,d in SD.G.nodes(data=True)) == 1
    
    
# TODO: add tests to verify the attractors are properly found
# TODO: add tests for a wider variety of networks