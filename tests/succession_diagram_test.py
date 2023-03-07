from nfvsmotifs.SuccessionDiagram import SuccessionDiagram
from biodivine_aeon import BooleanNetwork # type: ignore

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
    

def test_succession_diagram_structure2():
    bn = BooleanNetwork.from_bnet("""
        x1, x2
        x2, x1
        x3, x1 | !x4
        x4, x3
        x5, x2 & x6
        x6, x5
        """)

    SD = SuccessionDiagram(bn)
    assert SD.G.number_of_nodes() == 6 # NOTE: We are using a slightly different convention; 
                                       # we are merging if the fixed vars match AFTER percolation, 
                                       # whereas pystablemotifs merges BEFORE percolation
    assert SD.G.number_of_edges() == 7
    assert max(d['depth'] for n,d in SD.G.nodes(data=True)) == 2
    
# TODO: add tests to verify the attractors are properly found
# TODO: add tests for a wider variety of networks