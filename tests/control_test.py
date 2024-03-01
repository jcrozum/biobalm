from biodivine_aeon import BooleanNetwork

from balm.control import (
    Intervention,
    controls_are_equal,
    drivers_of_succession,
    succession_control,
    successions_to_target,
)
from balm.succession_diagram import SuccessionDiagram
from balm.types import BooleanSpace


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
    target_succession: list[BooleanSpace] = [
        {"S": 0},
        {"A": 0, "B": 0},
        {"C": 1, "D": 1},
    ]

    cs: list[list[BooleanSpace]] = [
        [{"S": 0}],
        [{"A": 0}, {"B": 0}],
        [{"C": 1}, {"D": 1}],
    ]

    drivers = drivers_of_succession(bn, target_succession)
    assert all([controls_are_equal(a, b) for a, b in zip(cs, drivers)])

    target_succession = [
        {"S": 0},
        {"C": 1, "D": 1},
        {"A": 0, "B": 0},
    ]

    cs = [[{"S": 0}], [{"C": 1}, {"D": 1}], [{"A": 0}, {"B": 0}]]

    drivers = drivers_of_succession(bn, target_succession)
    assert all([controls_are_equal(a, b) for a, b in zip(cs, drivers)])

    target_succession = [
        {"E": 1},
        {"S": 0},
        {"C": 1, "D": 1},
        {"A": 0, "B": 0},
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
            {"S": 0},
            {"A": 0, "B": 0},
            {"C": 1, "D": 1},
        ],
        [
            {"S": 0},
            {"C": 1, "D": 1},
            {"A": 0, "B": 0},
        ],
    ]
    target: BooleanSpace = {"S": 0, "E": 0, "A": 0, "B": 0, "C": 1, "D": 1}
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
    sd = SuccessionDiagram.from_rules(
        """
    S, S
    A, S | B
    B, A
    C, A | D
    D, C
    E, false
    """
    )
    target: BooleanSpace = {"S": 0, "E": 0, "A": 0, "B": 0, "C": 1, "D": 1}

    true_controls: list[list[list[BooleanSpace]]] = [
        [[{"S": 0}], [{"A": 0}, {"B": 0}], [{"C": 1}, {"D": 1}]],
        [[{"S": 0}], [{"C": 1}, {"D": 1}], [{"A": 0}, {"B": 0}]],
    ]

    true_successions: list[list[BooleanSpace]] = [
        [
            {"S": 0},
            {"A": 0, "B": 0},
            {"C": 1, "D": 1},
        ],
        [
            {"S": 0},
            {"C": 1, "D": 1},
            {"A": 0, "B": 0},
        ],
    ]

    true_interventions = [
        Intervention(c, "internal", s) for c, s in zip(true_controls, true_successions)
    ]

    interventions = succession_control(sd, target)

    assert len(interventions) == len(true_interventions)
    for intervention in interventions:
        assert intervention in true_interventions


def test_all_succession_control():
    sd = SuccessionDiagram.from_rules(
        """
    S, S
    A, S | B
    B, A
    C, A | D
    D, C
    E, false
    """
    )
    target: BooleanSpace = {"S": 0, "E": 0, "A": 0, "B": 0, "C": 1, "D": 1}

    true_controls: list[list[list[BooleanSpace]]] = [
        [[{"S": 0}], [{"A": 0}, {"B": 0}], [{"C": 1}, {"D": 1}]],
        [[{"S": 0}], [{"A": 1}, {"B": 1}, {"C": 1}, {"D": 1}], [{"A": 0}, {"B": 0}]],
    ]

    true_successions: list[list[BooleanSpace]] = [
        [
            {"S": 0},
            {"A": 0, "B": 0},
            {"C": 1, "D": 1},
        ],
        [
            {"S": 0},
            {"C": 1, "D": 1},
            {"A": 0, "B": 0},
        ],
    ]

    true_interventions = [
        Intervention(c, "all", s) for c, s in zip(true_controls, true_successions)
    ]

    interventions = succession_control(sd, target, strategy="all")

    assert len(interventions) == len(true_interventions)
    for intervention in interventions:
        assert intervention in true_interventions


def test_forbidden_drivers():
    sd = SuccessionDiagram.from_rules(
        """
    A, B & C
    B, A & C
    C, A & B
    """
    )
    target: BooleanSpace = {"A": 1, "B": 1, "C": 1}

    # Test with no forbidden drivers first
    true_controls: list[list[list[BooleanSpace]]] = [
        [[{"A": 1, "B": 1}, {"A": 1, "C": 1}, {"B": 1, "C": 1}]]
    ]
    true_successions: list[list[BooleanSpace]] = [[{"A": 1, "B": 1, "C": 1}]]

    true_interventions = [
        Intervention(c, "internal", s) for c, s in zip(true_controls, true_successions)
    ]

    interventions = succession_control(sd, target)

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

    interventions = succession_control(sd, target, forbidden_drivers=forbidden_drivers)

    assert len(interventions) == len(true_interventions)
    for intervention in interventions:
        assert intervention in true_interventions
        assert intervention.successful

    # Test with forbidden driver; case without solution
    forbidden_drivers: set[str] = set(("A", "B"))  # type: ignore

    true_controls = [[[]]]  # type: ignore
    true_successions = [[{"A": 1, "B": 1, "C": 1}]]

    true_interventions = [
        Intervention(c, "internal", s)
        for c, s in zip(true_controls, true_successions)  # type: ignore
    ]

    # do not show failed solution (default)
    interventions = succession_control(sd, target, forbidden_drivers=forbidden_drivers)
    assert len(interventions) == 0

    # show failed solution
    interventions = succession_control(
        sd, target, forbidden_drivers=forbidden_drivers, successful_only=False
    )

    assert len(interventions) == len(true_interventions)
    for intervention in interventions:
        assert intervention in true_interventions
        assert not intervention.successful


def test_size_restriction():
    sd = SuccessionDiagram.from_rules(
        """
    A, B & C
    B, A & C
    C, A & B
    """
    )
    target: BooleanSpace = {"A": 1, "B": 1, "C": 1}

    # Test with no restrictions
    true_controls: list[list[list[BooleanSpace]]] = [
        [[{"A": 1, "B": 1}, {"A": 1, "C": 1}, {"B": 1, "C": 1}]]
    ]
    true_successions: list[list[BooleanSpace]] = [[{"A": 1, "B": 1, "C": 1}]]

    true_interventions = [
        Intervention(c, "internal", s) for c, s in zip(true_controls, true_successions)
    ]

    interventions = succession_control(sd, target)

    assert len(interventions) == len(true_interventions)
    for intervention in interventions:
        assert intervention in true_interventions
        assert intervention.successful

    # Test with size restriction; no solution exists
    true_controls = [[[]]]  # type: ignore
    true_successions = [[{"A": 1, "B": 1, "C": 1}]]

    true_interventions = [
        Intervention(c, "internal", s)
        for c, s in zip(true_controls, true_successions)  # type: ignore
    ]

    # show the failed solution
    interventions = succession_control(
        sd, target, max_drivers_per_succession_node=1, successful_only=False
    )

    assert len(interventions) == len(true_interventions)
    for intervention in interventions:
        assert intervention in true_interventions
        assert not intervention.successful

    # do not show the failed solution (default)
    interventions = succession_control(
        sd, target, max_drivers_per_succession_node=1, successful_only=True
    )
    assert len(interventions) == 0
