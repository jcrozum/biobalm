from biodivine_aeon import BooleanNetwork

from nfvsmotifs.control import (
    Intervention,
    controls_are_equal,
    drivers_of_succession,
    succession_control,
    successions_to_target,
)
from nfvsmotifs.SuccessionDiagram import SuccessionDiagram


def test_intervention_equality_and_equivalence():
    i1 = Intervention([[{"Y": 0}, {"X": 0}]], "internal", [{}])
    i2 = Intervention([[{"X": 0}, {"Y": 0}]], "internal", [{}])
    assert i1 == i2

    i1 = Intervention([[{"Y": 0}, {"X": 0}]], "internal", [{"A": 0}])
    i2 = Intervention([[{"X": 0}, {"Y": 0}]], "internal", [{}])
    assert i1 != i2
    assert i1.is_equivalent(i2)

    i1 = Intervention([[{"Y": 0}, {"X": 0}]], "all", [{}])
    i2 = Intervention([[{"X": 0}, {"Y": 0}]], "all", [{}])
    assert i1 == i2

    i1 = Intervention([[{"Y": 0}, {"X": 0}]], "all", [{"A": 0}])
    i2 = Intervention([[{"X": 0}, {"Y": 0}]], "all", [{}])
    assert i1 != i2
    assert not i1.is_equivalent(i2)


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

    drivers = drivers_of_succession(bn, target_succession)
    assert all([controls_are_equal(a, b) for a, b in zip(cs, drivers)])

    target_succession = [
        {"S": 0, "E": 0},
        {"S": 0, "E": 0, "C": 1, "D": 1},
        {"S": 0, "E": 0, "A": 0, "B": 0, "C": 1, "D": 1},
    ]

    cs = [[{"S": 0}], [{"C": 1}, {"D": 1}], [{"A": 0}, {"B": 0}]]

    drivers = drivers_of_succession(bn, target_succession)
    assert all([controls_are_equal(a, b) for a, b in zip(cs, drivers)])

    target_succession = [
        {"E": 1},
        {"S": 0, "E": 1},
        {"S": 0, "E": 1, "C": 1, "D": 1},
        {"S": 0, "E": 1, "A": 0, "B": 0, "C": 1, "D": 1},
    ]

    cs = [[{"E": 1}], [{"S": 0}], [{"C": 1}, {"D": 1}], [{"A": 0}, {"B": 0}]]

    drivers = drivers_of_succession(bn, target_succession)
    assert all([controls_are_equal(a, b) for a, b in zip(cs, drivers)])


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
        tuple(frozenset(ts.items()) for ts in succession) for succession in successions
    )
    targets_hashed = set(
        tuple(frozenset(ts.items()) for ts in succession)
        for succession in target_successions
    )

    assert targets_hashed == successions_hashed


def test_internal_succession_control():
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
    target = {"S": 0, "E": 0, "A": 0, "B": 0, "C": 1, "D": 1}

    true_controls = [
        [[{"S": 0}], [{"A": 0}, {"B": 0}], [{"C": 1}, {"D": 1}]],
        [[{"S": 0}], [{"C": 1}, {"D": 1}], [{"A": 0}, {"B": 0}]],
    ]

    true_successions = [
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

    true_interventions = [
        Intervention(c, "internal", s) for c, s in zip(true_controls, true_successions)
    ]

    interventions = succession_control(bn, target)

    assert len(interventions) == len(true_interventions)
    for intervention in interventions:
        assert intervention in true_interventions


def test_all_succession_control():
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
    target = {"S": 0, "E": 0, "A": 0, "B": 0, "C": 1, "D": 1}

    true_controls = [
        [[{"S": 0}], [{"A": 0}, {"B": 0}], [{"C": 1}, {"D": 1}]],
        [[{"S": 0}], [{"A": 1}, {"B": 1}, {"C": 1}, {"D": 1}], [{"A": 0}, {"B": 0}]],
    ]

    true_successions = [
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

    true_interventions = [
        Intervention(c, "all", s) for c, s in zip(true_controls, true_successions)
    ]

    interventions = succession_control(bn, target, strategy="all")

    assert len(interventions) == len(true_interventions)
    for intervention in interventions:
        assert intervention in true_interventions


def test_forbidden_drivers():
    bn = BooleanNetwork.from_bnet(
        """
    A, B & C
    B, A & C
    C, A & B
    """
    )
    target = {"A": 1, "B": 1, "C": 1}

    # Test with no forbidden drivers first
    true_controls = [[[{"A": 1, "B": 1}, {"A": 1, "C": 1}, {"B": 1, "C": 1}]]]
    true_successions = [[{"A": 1, "B": 1, "C": 1}]]

    true_interventions = [
        Intervention(c, "internal", s) for c, s in zip(true_controls, true_successions)
    ]

    interventions = succession_control(bn, target)

    assert len(interventions) == len(true_interventions)
    for intervention in interventions:
        assert intervention in true_interventions
        assert intervention.successful

    # Test with forbidden driver; case with solution
    forbidden_drivers = set("A")

    true_controls = [[[{"B": 1, "C": 1}]]]
    true_successions = [[{"A": 1, "B": 1, "C": 1}]]

    true_interventions = [
        Intervention(c, "internal", s) for c, s in zip(true_controls, true_successions)
    ]

    interventions = succession_control(bn, target, forbidden_drivers=forbidden_drivers)

    assert len(interventions) == len(true_interventions)
    for intervention in interventions:
        assert intervention in true_interventions
        assert intervention.successful

    # Test with forbidden driver; case without solution
    forbidden_drivers = set(("A", "B"))

    true_controls = [[[]]]  # type: ignore
    true_successions = [[{"A": 1, "B": 1, "C": 1}]]

    true_interventions = [
        Intervention(c, "internal", s) for c, s in zip(true_controls, true_successions)  # type: ignore
    ]

    interventions = succession_control(bn, target, forbidden_drivers=forbidden_drivers)

    assert len(interventions) == len(true_interventions)
    for intervention in interventions:
        assert intervention in true_interventions
        assert not intervention.successful


def test_size_restriction():
    bn = BooleanNetwork.from_bnet(
        """
    A, B & C
    B, A & C
    C, A & B
    """
    )
    target = {"A": 1, "B": 1, "C": 1}

    # Test with no restrictions
    true_controls = [[[{"A": 1, "B": 1}, {"A": 1, "C": 1}, {"B": 1, "C": 1}]]]
    true_successions = [[{"A": 1, "B": 1, "C": 1}]]

    true_interventions = [
        Intervention(c, "internal", s) for c, s in zip(true_controls, true_successions)
    ]

    interventions = succession_control(bn, target)

    assert len(interventions) == len(true_interventions)
    for intervention in interventions:
        assert intervention in true_interventions
        assert intervention.successful

    # Test with size restriction; no solution exists
    true_controls = [[[]]]  # type: ignore
    true_successions = [[{"A": 1, "B": 1, "C": 1}]]

    true_interventions = [
        Intervention(c, "internal", s) for c, s in zip(true_controls, true_successions)  # type: ignore
    ]

    interventions = succession_control(bn, target, max_drivers_per_succession_node=1)

    assert len(interventions) == len(true_interventions)
    for intervention in interventions:
        assert intervention in true_interventions
        assert not intervention.successful
