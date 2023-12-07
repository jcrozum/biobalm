from biodivine_aeon import BooleanNetwork  # type:ignore

from balm.motif_avoidant import _filter_candidates  # type: ignore
from balm.motif_avoidant import _Pint_reachability  # type: ignore
from balm.motif_avoidant import _preprocess_candidates  # type: ignore
from balm.petri_net_translation import network_to_petrinet
from balm.state_utils import state_list_to_bdd


def test_preprocessing_ssf_not_optimal():
    bn = BooleanNetwork.from_bnet(
        """
        x1, (x1 & x2) | (!x1 & !x2)
        x2, (x1 & x2) | (!x1 & !x2)
    """
    )

    s0 = {"x1": 0, "x2": 0}
    s1 = {"x1": 0, "x2": 1}
    s2 = {"x1": 1, "x2": 0}
    # s3 = {"x1": 1, "x2": 1}

    """
        This BN has one minimal trap space: 11.
        It has one motif-avoidant attractor: {00, 01, 10}.
        It has one time-reversal minimal trap space: 11, which is not self-negating.
        Assume that terminal_restriction_space = {0, 1}^2 - {11} = {00, 01, 10}.
    """

    terminal_restriction_space = state_list_to_bdd([s0, s1, s2])

    """
        The minimum NFVS is {x1, x2}.
        If b_1 = 0 and b_2 = 0, then F = [00].
        If b_1 = 1 and b_2 = 1, then F = [01, 10].
    """

    # candidates_F = {00}
    candidates_F = [s0]
    candidates_F = _preprocess_candidates(
        bn, candidates_F, terminal_restriction_space, 1000
    )
    assert len(candidates_F) == 1

    # candidates_F = {01, 10}
    candidates_F = [s1, s2]
    candidates_F = _preprocess_candidates(
        bn, candidates_F, terminal_restriction_space, 1000
    )
    assert len(candidates_F) == 1


def test_preprocessing_ssf_optimal():
    bn = BooleanNetwork.from_bnet(
        """
        A, !B
        B, !A
        C, A | B
    """
    )

    s0 = {"A": 0, "B": 0, "C": 0}
    # s1 = {"A": 0, "B": 0, "C": 1}
    s2 = {"A": 0, "B": 1, "C": 0}
    # s3 = {"A": 0, "B": 1, "C": 1}
    s4 = {"A": 1, "B": 0, "C": 0}
    # s5 = {"A": 1, "B": 0, "C": 1}
    # s6 = {"A": 1, "B": 1, "C": 0}
    s7 = {"A": 1, "B": 1, "C": 1}

    """
        This BN has two minimal trap spaces: 101 + 011.
        It has no motif-avoidant attractor.
        It has two time-reversal minimal trap spaces: 110 + 001. Both are self-negating.
        Assume that terminal_restriction_space = {0, 1}^3 - {101, 011, 110, 001} = {000, 010, 100, 111}.
    """
    terminal_restriction_space = state_list_to_bdd([s0, s2, s4, s7])

    """
        The minimum NFVS is {A, B}.
        Assume that b_A = 0 and b_B = 0.
        Then F = [000].
    """

    candidates_F = [s0]
    candidates_F = _preprocess_candidates(
        bn, candidates_F, terminal_restriction_space, 1000
    )
    assert len(candidates_F) == 0


def test_ABNReach_current_version():
    bn = BooleanNetwork.from_bnet(
        """
        x1, (x1 & x2) | (!x1 & !x2)
        x2, (x1 & x2) | (!x1 & !x2)
        x3, x3 | !x3
    """
    )

    s0 = {"x1": 0, "x2": 0, "x3": 1}
    s1 = {"x1": 0, "x2": 1, "x3": 1}
    # s2 = {"x1": 1, "x2": 0, "x3": 1}
    s3 = {"x1": 1, "x2": 1, "x3": 1}

    petri_net = network_to_petrinet(bn)

    joint_target_set = state_list_to_bdd([s3])
    is_reachable = _Pint_reachability(petri_net, s0, joint_target_set)
    assert (
        is_reachable is False
    )  # 00 does not reach 11, Pint cannot determinem but Mole can

    joint_target_set = state_list_to_bdd([s0])
    is_reachable = _Pint_reachability(petri_net, s3, joint_target_set)
    assert is_reachable is False  # 11 does not reach 00

    joint_target_set = state_list_to_bdd([s1])
    is_reachable = _Pint_reachability(petri_net, s0, joint_target_set)
    assert is_reachable is True  # 00 reaches 01

    joint_target_set = state_list_to_bdd([s1, s3])
    is_reachable = _Pint_reachability(petri_net, s0, joint_target_set)
    assert is_reachable is True  # 00 reaches 01


def test_FilteringProcess():
    bn = BooleanNetwork.from_bnet(
        """
        x1, (x1 & x2) | (!x1 & !x2)
        x2, (x1 & x2) | (!x1 & !x2)
    """
    )

    s0 = {"x1": 0, "x2": 0}
    s1 = {"x1": 0, "x2": 1}
    s2 = {"x1": 1, "x2": 0}
    # s3 = {"x1": 1, "x2": 1}

    terminal_res_space = state_list_to_bdd([s0, s1, s2])
    petri_net = network_to_petrinet(bn)

    F = [s1, s2]  # Candidate set after finishing preprocessing
    list_motif_avoidant_atts = _filter_candidates(petri_net, F, terminal_res_space)
    assert len(list_motif_avoidant_atts) == 1  # a motif-avoidant attractor {00, 01, 10}
