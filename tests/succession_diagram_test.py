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

def test_succession_diagram_structure3():
    bn = BooleanNetwork.from_bnet("""
        AUXINS,  AUXINS | !AUXINS
        SHR,     SHR
        ARF,     !IAA
        IAA,     !AUXINS
        JKD,     SHR & SCR
        MGP,     !WOX & SHR & SCR
        SCR,     SHR & SCR & !MGP | SHR & SCR & JKD
        WOX,     WOX & SHR & SCR & ARF | SHR & SCR & !MGP & ARF
        PLT,     ARF
        """)
    
    atts_true = [
        {'ARF': 1, 'AUXINS': 1, 'IAA': 0, 'JKD': 0, 'MGP': 0, 'PLT': 1, 'SCR': 0, 'SHR': 0, 'WOX': 0},
        {'ARF': 1, 'AUXINS': 1, 'IAA': 0, 'JKD': 0, 'MGP': 0, 'PLT': 1, 'SCR': 0, 'SHR': 1, 'WOX': 0},
        {'ARF': 1, 'AUXINS': 1, 'IAA': 0, 'JKD': 1, 'MGP': 0, 'PLT': 1, 'SCR': 1, 'SHR': 1, 'WOX': 1},
        {'ARF': 1, 'AUXINS': 1, 'IAA': 0, 'JKD': 1, 'MGP': 1, 'PLT': 1, 'SCR': 1, 'SHR': 1, 'WOX': 0},
    ]

    SD = SuccessionDiagram(bn)
    
    atts_found = [d['attractors'][0] for n,d in SD.G.nodes(data=True) if len(d['attractors']) > 0]
    
    assert len(atts_true) == len(atts_found)
    assert atts_found[0] in atts_true
    assert atts_found[1] in atts_true
    assert atts_found[2] in atts_true
    assert atts_found[3] in atts_true