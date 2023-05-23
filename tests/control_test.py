from biodivine_aeon import BooleanNetwork

from nfvsmotifs.control import drivers_of_succession, successions_to_target
from nfvsmotifs.SuccessionDiagram import SuccessionDiagram


def test_basic_succession_control():
    bn = BooleanNetwork.from_bnet(
        """
    S, S
    A, S | B
    B, A
    C, A | D
    D, C
    E, false
    """
    )
    target_succession = [
        {"S": 0, "E": 0},
        {"S": 0, "E": 0, "A": 0, "B": 0},
        {"S": 0, "E": 0, "A": 0, "B": 0, "C": 1, "D": 1},
    ]

    cs = [[{"S": 0}], [{"A": 0}, {"B": 0}], [{"C": 1}, {"D": 1}]]

    # NOTE: this test is not perfect right now because
    # the ordering of the inner lists should not matter
    assert drivers_of_succession(bn, target_succession) == cs

    target_succession = [
        {"S": 0, "E": 0},
        {"S": 0, "E": 0, "C": 1, "D": 1},
        {"S": 0, "E": 0, "A": 0, "B": 0, "C": 1, "D": 1},
    ]

    cs = [[{"S": 0}], [{"C": 1}, {"D": 1}], [{"A": 0}, {"B": 0}]]

    # NOTE: this test is not perfect right now because
    # the ordering of the inner lists should not matter
    assert drivers_of_succession(bn, target_succession) == cs

    target_succession = [
        {"E": 1},
        {"S": 0, "E": 1},
        {"S": 0, "E": 1, "C": 1, "D": 1},
        {"S": 0, "E": 1, "A": 0, "B": 0, "C": 1, "D": 1},
    ]

    cs = [[{"E": 1}], [{"S": 0}], [{"C": 1}, {"D": 1}], [{"A": 0}, {"B": 0}]]

    # NOTE: this test is not perfect right now because
    # the ordering of the inner lists should not matter
    assert drivers_of_succession(bn, target_succession) == cs


def test_basic_succession_finding():
    bn = BooleanNetwork.from_bnet(
        """
    S, S
    A, S | B
    B, A
    C, A | D
    D, C
    E, false
    """
    )
    target_successions = [
        [
            {"S": 0, "E": 0},
            {"S": 0, "E": 0, "A": 0, "B": 0},
            {"S": 0, "E": 0, "A": 0, "B": 0, "C": 1, "D": 1},
        ],
        [
            {"S": 0, "E": 0},
            {"S": 0, "E": 0, "C": 1, "D": 1},
            {"S": 0, "E": 0, "A": 0, "B": 0, "C": 1, "D": 1},
        ],
    ]
    target = {"S": 0, "E": 0, "A": 0, "B": 0, "C": 1, "D": 1}
    succession_diagram = SuccessionDiagram(bn)

    successions = successions_to_target(succession_diagram, target)

    successions_hashed = set(
        frozenset(ts.items()) for succession in successions for ts in succession
    )
    targets_hashed = set(
        frozenset(ts.items()) for succession in target_successions for ts in succession
    )

    assert targets_hashed == successions_hashed
