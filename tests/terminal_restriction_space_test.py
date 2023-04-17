from biodivine_aeon import BooleanNetwork # type: ignore
from nfvsmotifs.trappist_core import trappist
from nfvsmotifs.terminal_restriction_space import *

def test_tr_trap_spaces():
    """
    TODO: need to make a test using real models
    """

    bn = BooleanNetwork.from_bnet("""
        A, B
        B, A | B
    """)

    tr_trap_spaces = trappist(bn, problem="max", reverse_time=True)

    assert {'A': 0, 'B': 1} in tr_trap_spaces
    assert {'B': 0} in tr_trap_spaces
    assert [{'A': 0, 'B': 1}] == get_self_neg_tr_trap_spaces(bn)

def test_get_terminal_restriction_space():
    network = BooleanNetwork.from_bnet("""
    A, !A & !B | C
    B, !A & !B | C
    C, A & B
    """)
    stable_motifs = [{'A':1, 'B':1, 'C':1}]

    TRS = get_terminal_restriction_space(stable_motifs, network, use_single_node_drivers = True, use_tr_trapspaces = True)

    assert TRS == state_list_to_bdd([{'A':0,'B':0,'C':0},{'A':1,'B':0,'C':0},{'A':0,'B':1,'C':0}])