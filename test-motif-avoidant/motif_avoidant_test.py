from pyeda.boolalg import boolfunc # type:ignore
from pyeda.boolalg.bdd import bddvar, expr2bdd, BinaryDecisionDiagram # type:ignore
from pyeda.boolalg.expr import expr # type:ignore

from typing import List, Set, Dict # type: ignore

from biodivine_aeon import BooleanNetwork # type:ignore
from nfvsmotifs.motif_avoidant import PreprocessingSSF, FilteringProcess, ABNReach_current_version # type:ignore
from nfvsmotifs.pyeda_utils import aeon_to_pyeda # type:ignore
from nfvsmotifs.state_utils import state_2_bdd, list_state_2_bdd, eval_function, is_member_bdd # type:ignore
from nfvsmotifs.petri_net_translation import network_to_petrinet # type:ignore

def test_preprocessing_ssf_not_optimal():
    bn = BooleanNetwork.from_bnet("""
        x1, (x1 & x2) | (!x1 & !x2)
        x2, (x1 & x2) | (!x1 & !x2)
    """)
    
    s0 = {'x1': 0, 'x2': 0}
    s1 = {'x1': 0, 'x2': 1}
    s2 = {'x1': 1, 'x2': 0}
    s3 = {'x1': 1, 'x2': 1}

    """
        This ABN has one minimal trap spaces: 11.
        It has one motif-avoidant attractor: {00, 01, 10}.
        It has one time-reversal minimal trap space: 11, which is not self-negative.
        Assume that terminal_res_space = {0, 1}^2 - {11} = {00, 01, 10}.
    """

    terminal_res_space = list_state_2_bdd([s0, s1, s2])

    """
        The minimum NFVS is {x1, x2}.
        If b_1 = 0 and b_2 = 0, then F = [00].
        If b_1 = 1 and b_2 = 1, then F = [01, 10].
    """

    # F = {00}
    F = [s0]
    F = PreprocessingSSF(bn, F, terminal_res_space, I_MAX=1000)
    assert len(F) == 1

    # F = {01, 10}
    F = [s1, s2]
    F = PreprocessingSSF(bn, F, terminal_res_space, I_MAX=1000)
    assert len(F) == 1
    

def test_preprocessing_ssf_optimal():
    bn = BooleanNetwork.from_bnet("""
        A, !B
        B, !A
        C, A | B
    """)
    
    s0 = {'A': 0, 'B': 0, 'C': 0}
    s1 = {'A': 0, 'B': 0, 'C': 1}
    s2 = {'A': 0, 'B': 1, 'C': 0}
    s3 = {'A': 0, 'B': 1, 'C': 1}
    s4 = {'A': 1, 'B': 0, 'C': 0}
    s5 = {'A': 1, 'B': 0, 'C': 1}
    s6 = {'A': 1, 'B': 1, 'C': 0}
    s7 = {'A': 1, 'B': 1, 'C': 1}

    """
        This ABN has two minimal trap spaces: 101 + 011.
        It has no motif-avoidant attractor.
        It has two time-reversal minimal trap spaces: 110 + 001. Both are self-negative ones.
        Assume that terminal_res_space = {0, 1}^3 - {101, 011, 110, 001} = {000, 010, 100, 111}.
    """
    terminal_res_space = list_state_2_bdd([s0, s2, s4, s7])

    """
        The minimum NFVS is {A, B}.
        Assume that b_A = 0 and b_B = 0.
        Then F = [000].
    """

    F = [s0]
    F = PreprocessingSSF(bn, F, terminal_res_space, I_MAX=1000)
    assert len(F) == 0


def test_ABNReach_current_version():
    bn = BooleanNetwork.from_bnet("""
        x1, (x1 & x2) | (!x1 & !x2)
        x2, (x1 & x2) | (!x1 & !x2)
    """)
    
    s0 = {'x1': 0, 'x2': 0}
    s1 = {'x1': 0, 'x2': 1}
    s2 = {'x1': 1, 'x2': 0}
    s3 = {'x1': 1, 'x2': 1}

    petri_net = network_to_petrinet(bn)

    joint_target_set = list_state_2_bdd([s3])
    is_reachable = ABNReach_current_version(bn, petri_net, s0, joint_target_set)
    assert is_reachable ==  False # 00 does not reach 11, Pint cannot determinem but Mole can

    joint_target_set = list_state_2_bdd([s0])
    is_reachable = ABNReach_current_version(bn, petri_net, s3, joint_target_set)
    assert is_reachable == False # 11 does not reach 00

    joint_target_set = list_state_2_bdd([s1])
    is_reachable = ABNReach_current_version(bn, petri_net, s0, joint_target_set)
    assert is_reachable == True # 00 reaches 01

    joint_target_set = list_state_2_bdd([s1, s3])
    is_reachable = ABNReach_current_version(bn, petri_net, s0, joint_target_set)
    assert is_reachable == True # 00 reaches 01


def test_FilteringProcess():
    bn = BooleanNetwork.from_bnet("""
        x1, (x1 & x2) | (!x1 & !x2)
        x2, (x1 & x2) | (!x1 & !x2)
    """)
    
    s0 = {'x1': 0, 'x2': 0}
    s1 = {'x1': 0, 'x2': 1}
    s2 = {'x1': 1, 'x2': 0}
    s3 = {'x1': 1, 'x2': 1}

    terminal_res_space = list_state_2_bdd([s0, s1, s2])
    petri_net = network_to_petrinet(bn)

    F = [s1, s2] # Candidate set after finishing PreprocessingSSF
    list_motif_avoidant_atts = FilteringProcess(bn, petri_net, F, terminal_res_space)
    assert len(list_motif_avoidant_atts) == 1 # a motif-avoidant attractor {00, 01, 10}


