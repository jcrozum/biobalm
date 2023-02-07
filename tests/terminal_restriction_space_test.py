from biodivine_aeon import BooleanNetwork   
from nfvsmotifs.trappist_core import trappist
from nfvsmotifs.terminal_restriction_space import *
from nfvsmotifs.aeon_utils import remove_static_constraints

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
