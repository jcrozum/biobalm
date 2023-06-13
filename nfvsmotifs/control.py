from __future__ import annotations

from itertools import combinations, product
from typing import cast

import networkx as nx
from biodivine_aeon import BooleanNetwork

from nfvsmotifs.space_utils import percolate_space
from nfvsmotifs.SuccessionDiagram import SuccessionDiagram

SuccessionType = list[dict[str, int]]
ControlType = list[dict[str, int]]


def succession_control(
    bn: BooleanNetwork,
    target: dict[str, int],
    strategy: str = "internal",
    succession_diagram: SuccessionDiagram | None = None,
) -> list[tuple[list[ControlType], SuccessionType]]:
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

    Returns
    -------
    list[tuple[ControlType, SuccessionType]]
        A list of control strategies. Each control strategy is represented as
        tuple containing a list of state dictionaries and the succession from
        which it was derived.
    """
    interventions: list[tuple[list[ControlType], SuccessionType]] = []

    if succession_diagram is None:
        succession_diagram = SuccessionDiagram(bn)

    successions = successions_to_target(
        succession_diagram, target=target, expand_diagram=True
    )

    for succession in successions:
        controls = drivers_of_succession(bn, succession, strategy=strategy)
        interventions.append((controls, succession))
        # for cs in controls:
        #     interventions.append((cs, succession))

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
            find_drivers(bn, ts, strategy=strategy, assume_fixed=assume_fixed)
        )
        ldoi, _ = percolate_space(bn, ts | assume_fixed, strict_percolation=False)
        assume_fixed.update(ldoi)

    return control_strategies


def find_drivers(
    bn: BooleanNetwork,
    target_trap_space: dict[str, int],
    strategy: str = "internal",
    assume_fixed: dict[str, int] | None = None,
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

    Returns
    -------
    ControlType
        A list of internal driver sets, represented as state dictionaries.
    """
    if assume_fixed is None:
        assume_fixed = {}

    target_trap_space_inner = {
        k: v for k, v in target_trap_space.items() if k not in assume_fixed
    }

    if strategy == "internal":
        driver_pool = list(target_trap_space_inner)
    elif strategy == "all":
        driver_pool = [bn.get_variable_name(id) for id in bn.variables()]
    else:
        raise ValueError("Unknown driver search strategy")

    max_drivers = len(target_trap_space_inner)
    drivers: ControlType = []
    for driver_set_size in range(1, max_drivers + 1):
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
