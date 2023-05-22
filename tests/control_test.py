from biodivine_aeon import BooleanNetwork

from nfvsmotifs.control import drivers_of_succession


def test_basic_succession_control():
    bn = BooleanNetwork.from_bnet(
        """
    S, S
    A, S | B
    B, A
    C, A | D
    D, C
    """
    )
    target_succession = [
        {"S": 0},
        {"S": 0, "A": 0, "B": 0},
        {"S": 0, "A": 0, "B": 0, "C": 1, "D": 1},
    ]

    cs = [[{"S": 0}], [{"A": 0}, {"B": 0}], [{"C": 1}, {"D": 1}]]

    # NOTE: this test is not perfect right now because the ordering of the inner lists should not matter
    assert drivers_of_succession(bn, target_succession) == cs
