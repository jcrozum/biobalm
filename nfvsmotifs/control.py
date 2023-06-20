from __future__ import annotations

from itertools import combinations, product
from typing import cast

import networkx as nx
from biodivine_aeon import BooleanNetwork

from nfvsmotifs.space_utils import percolate_space
from nfvsmotifs.SuccessionDiagram import SuccessionDiagram

SuccessionType = list[dict[str, int]]  # sequence of stable motifs
ControlType = list[dict[str, int]]  # ways of locking in an individual stable motif


def controls_are_equal(a: ControlType, b: ControlType) -> bool:
    return set(frozenset(x.items()) for x in a) == set(frozenset(x.items()) for x in b)


class Intervention:
    def __init__(
        self, control: list[ControlType], strategy: str, succession: SuccessionType
    ):
        self._control = control
        self._strategy = strategy
        self._succession = succession
        self._successful = not any(not c for c in control)

    @property
    def control(self):
        return self._control

    @property
    def strategy(self):
        return self._strategy

    @property
    def succession(self):
        return self._succession

    @property
    def successful(self):
        return self._successful

    def is_equivalent(self, other: Intervention) -> bool:
        if self.strategy != other.strategy:
            return False

        # if using external drivers, the succession matters because it
        # determines how long you have to maintain temporary controls
        if self.strategy == "all":
            if self.succession != other.succession:
                return False

        if len(self.control) != len(other.control):
            return False

        for d1, d2 in zip(self.control, other.control):
            if not controls_are_equal(d1, d2):
                return False

        return True

    def __eq__(self, other: object):
        if not isinstance(other, Intervention):
            return False

        # if the strategy is "all", then is_equivalent will handle the
        # succession comparison
        if self.strategy != "all":
            if self.succession != other.succession:
                return False

        if not self.is_equivalent(other):
            return False

        return True

    def __repr__(self):
        return f"Intervention({self.control}, {self.strategy}, {self.succession})"

    def __str__(self):
        succession_string = (
            "operating on\n" + "\n".join(map(str, self.succession)) + "\noverride\n"
        )
        if self.strategy == "internal":
            return succession_string + " and \n".join(
                f"({' or '.join(map(str,motif_control))})"
                for motif_control in self.control
            )
        elif self.strategy == "all":
            return succession_string + "temporarily, and then \n".join(
                f"({' or '.join(map(str,motif_control))})"
                for motif_control in self.control
            )
        else:
            return "unknown strategy: " + self.__repr__()


def succession_control(
    bn: BooleanNetwork,
    target: dict[str, int],
    strategy: str = "internal",
    succession_diagram: SuccessionDiagram | None = None,
    max_drivers_per_succession_node: int | None = None,
    forbidden_drivers: set[str] | None = None,
) -> list[Intervention]:
    """_summary_

    Parameters
    ----------
    bn : BooleanNetwork
        The network to analyze, which contains the Boolean update functions.
    target : dict[str, int]
        The target subspace.
    strategy : str, optional
        The searching strategy to use to look for driver nodes. Options are
        'internal' (default), 'all'.
    succession_diagram : SuccessionDiagram | None, optional
        The succession diagram from which successions will be extracted. If
        `None`, then a succession diagram will be generated from `bn`.
    max_drivers_per_succession_node: int | None = None,
        The maximum number of drivers that will be tested for a succession
        diagram node. If `None`, then a number of drivers up to the size of the
        succession diagram node's stable motif will be tested
    forbidden_drivers: set[str] | None
        A set of forbidden drivers that will not be overridden for control. If
        `None`, then all nodes are candidates for control.

    Returns
    -------
    list[Intervention]
        A list of control intervention objects.
    """
    interventions: list[Intervention] = []

    if succession_diagram is None:
        succession_diagram = SuccessionDiagram(bn)

    successions = successions_to_target(
        succession_diagram, target=target, expand_diagram=True
    )

    for succession in successions:
        controls = drivers_of_succession(
            bn,
            succession,
            strategy=strategy,
            max_drivers_per_succession_node=max_drivers_per_succession_node,
            forbidden_drivers=forbidden_drivers,
        )
        interventions.append(Intervention(controls, strategy, succession))

    return interventions


def successions_to_target(
    succession_diagram: SuccessionDiagram,
    target: dict[str, int],
    expand_diagram: bool = True,
) -> list[SuccessionType]:
    """Find lists of nested trap spaces (successions) that lead to the
    specified target subspace.

    Parameters
    ----------
    succession_diagram : SuccessionDiagram
        The succession diagram from which successions will be extracted.
    target : dict[str, int]
        The target subspace.
    expand_diagram: bool
        Whether to ensure that the succession diagram is expanded enough to
        capture all paths to the target (default: True).

    Returns
    -------
    list[SuccessionType]
        A list of successions, where each succession is a list of sequentially
        nested trap spaces that specify the target.
    """
    successions: list[SuccessionType] = []

    # expand the succession_diagram toward the target
    if expand_diagram:
        succession_diagram.expand_node(
            succession_diagram.root(),
            depth_limit=None,
            node_limit=None,
            to_target=target,
        )

    for s in cast(list[int], succession_diagram.G.nodes()):
        fixed_vars = cast(dict[str, int], succession_diagram.G.nodes[s]["fixed_vars"])
        is_consistent = not any(
            k in target and target[k] != v for k, v in fixed_vars.items()
        )
        is_last_needed = set(target) <= set(fixed_vars)

        if not is_consistent or not is_last_needed:
            continue

        for path in cast(
            list[list[int]],
            nx.all_simple_paths(  # type: ignore
                succession_diagram.G,
                source=succession_diagram.root(),
                target=s,
            ),
        ):
            succession = [
                cast(dict[str, int], succession_diagram.G.edges[x, y]["motif"])
                for x, y in zip(path[:-1], path[1:])
            ]
            successions.append(succession)

    return successions


def drivers_of_succession(
    bn: BooleanNetwork,
    succession: list[dict[str, int]],
    strategy: str = "internal",
    max_drivers_per_succession_node: int | None = None,
    forbidden_drivers: set[str] | None = None,
) -> list[ControlType]:
    """Find driver nodes of a list of sequentially nested trap spaces

    Parameters
    ----------
    bn : BooleanNetwork
        The network to analyze, which contains the Boolean update functions.
    succession : list[dict[str, int]]
        A list of sequentially nested trap spaces that specify the target.
    strategy: str
        The searching strategy to use to look for driver nodes. Options are
        'internal' (default), 'all'.
    max_drivers_per_succession_node: int | None = None,
        The maximum number of drivers that will be tested for a succession
        diagram node. If `None`, then a number of drivers up to the size of the
        succession diagram node's stable motif will be tested
    forbidden_drivers: set[str] | None
        A set of forbidden drivers that will not be overridden for control. If
        `None`, then all nodes are candidates for control.

    Returns
    -------
    list[ControlType]
        A list of controls. Each control is a list of lists of driver sets,
        represented as state dictionaries. Each list item corresponds to a list
        of drivers for the corresponding trap space in the succession.
    """
    control_strategies: list[ControlType] = []
    assume_fixed: dict[str, int] = {}
    for ts in succession:
        control_strategies.append(
            find_drivers(
                bn,
                ts,
                strategy=strategy,
                assume_fixed=assume_fixed,
                max_drivers_per_succession_node=max_drivers_per_succession_node,
                forbidden_drivers=forbidden_drivers,
            )
        )
        ldoi, _ = percolate_space(bn, ts | assume_fixed, strict_percolation=False)
        assume_fixed.update(ldoi)

    return control_strategies


def find_drivers(
    bn: BooleanNetwork,
    target_trap_space: dict[str, int],
    strategy: str = "internal",
    assume_fixed: dict[str, int] | None = None,
    max_drivers_per_succession_node: int | None = None,
    forbidden_drivers: set[str] | None = None,
) -> ControlType:
    """Finds drives of a given target trap space

    Parameters
    ----------
    bn : BooleanNetwork
        The network to analyze, which contains the Boolean update functions.
    target_trap_space : dict[str, int]
        The trap space we want to find drivers for.
    strategy: str
        The searching strategy to use to look for driver nodes. Options are
        'internal' (default), 'all'.
    assume_fixed: dict[str,int] | None
        A dictionary of fixed variables that should be assumed to be fixed.
    max_drivers_per_succession_node: int | None = None,
        The maximum number of drivers that will be tested for a succession
        diagram node. If `None`, then a number of drivers up to the size of the
        succession diagram node's stable motif will be tested
    forbidden_drivers: set[str] | None
        A set of forbidden drivers that will not be overridden for control. If
        `None`, then all nodes are candidates for control.

    Returns
    -------
    ControlType
        A list of internal driver sets, represented as state dictionaries. If
        empty, then no drivers are found. This can happen if
        `max_drivers_per_succession_node` is not `None`, or if all controls
        require nodes in `forbidden_drivers`.
    """
    if assume_fixed is None:
        assume_fixed = {}
    if forbidden_drivers is None:
        forbidden_drivers = set()

    target_trap_space_inner = {
        k: v for k, v in target_trap_space.items() if k not in assume_fixed
    }

    if strategy == "internal":
        driver_pool = set(target_trap_space_inner) - forbidden_drivers
    elif strategy == "all":
        driver_pool = (
            set(bn.get_variable_name(id) for id in bn.variables()) - forbidden_drivers
        )
    else:
        raise ValueError("Unknown driver search strategy")

    if max_drivers_per_succession_node is None:
        max_drivers_per_succession_node = len(target_trap_space_inner)

    drivers: ControlType = []
    for driver_set_size in range(max_drivers_per_succession_node + 1):
        for driver_set in combinations(driver_pool, driver_set_size):
            if any(set(d) <= set(driver_set) for d in drivers):
                continue

            if strategy == "internal":
                driver_dict = {k: target_trap_space_inner[k] for k in driver_set}
                ldoi, _ = percolate_space(
                    bn, driver_dict | assume_fixed, strict_percolation=False
                )
                if target_trap_space.items() <= ldoi.items():
                    drivers.append(driver_dict)
            elif strategy == "all":
                for vals in product([0, 1], repeat=driver_set_size):
                    driver_dict = {
                        driver: value for driver, value in zip(driver_set, vals)
                    }
                    ldoi, _ = percolate_space(
                        bn, driver_dict | assume_fixed, strict_percolation=False
                    )
                    if target_trap_space.items() <= ldoi.items():
                        drivers.append(driver_dict)
    return drivers
