from biodivine_aeon import BooleanNetwork  # type: ignore

from nfvsmotifs.terminal_restriction_space import (
    get_self_neg_tr_trap_spaces,
    get_terminal_restriction_space,
    state_list_to_bdd,
)
from nfvsmotifs.trappist_core import trappist


def test_tr_trap_spaces():
    """
    TODO: need to make a test using real models
    """

    bn = BooleanNetwork.from_bnet(
        """
        A, B
        B, A | B
    """
    )

    tr_trap_spaces = trappist(bn, problem="max", reverse_time=True)

    assert {"A": 0, "B": 1} in tr_trap_spaces
    assert {"B": 0} in tr_trap_spaces
    assert [{"A": 0, "B": 1}] == get_self_neg_tr_trap_spaces(bn)


def test_get_terminal_restriction_space():
    network = BooleanNetwork.from_bnet(
        """
    A, !A & !B | C
    B, !A & !B | C
    C, A & B
    """
    )
    stable_motifs = [{"A": 1, "B": 1, "C": 1}]

    trs = get_terminal_restriction_space(
        stable_motifs,
        network,
        ensure_subspace={},
        use_single_node_drivers=True,
        use_tr_trapspaces=True,
    )

    assert trs == state_list_to_bdd(
        [{"A": 0, "B": 0, "C": 0}, {"A": 1, "B": 0, "C": 0}, {"A": 0, "B": 1, "C": 0}]
    )


def test_get_terminal_restriction_space2():
    network = BooleanNetwork.from_bnet(
        """
    A, !D | (A & !B & C)
    B, E & !(A & !B & C)
    C, F | (A & !B & C)
    D, C
    E, A
    F, B
    """
    )
    stable_motifs = [{"A": 1, "B": 0, "C": 1}]

    trs = get_terminal_restriction_space(
        stable_motifs,
        network,
        ensure_subspace={},
        use_single_node_drivers=True,
        use_tr_trapspaces=True,
    )

    assert trs == state_list_to_bdd(
        [
            {"A": 0},
            {"B": 1},
            {"C": 0},
        ]
    )
