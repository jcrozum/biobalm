from __future__ import annotations

from itertools import combinations

from biodivine_aeon import BooleanNetwork

from nfvsmotifs.space_utils import percolate_space


def drivers_of_succession(
    bn: BooleanNetwork,
    succession: list[dict[str, int]],
) -> list[list[dict[str, int]]]:
    """Find driver nodes of a list of sequentially nested trap spaces

    Parameters
    ----------
    bn : BooleanNetwork
        The network to analyze, which contains the Boolean update functions.
    succession : list[dict[str, int]]
        A list of sequentially nested trap spaces that specify the target.

    Returns
    -------
    list[list[dict[str, int]]]
        A list of lists of internal driver sets, represented as state dictionaries.
        Each list item corresponds to a list of drivers for the corresponding trap
        space in the succession.
    """
    control_strategies: list[list[dict[str, int]]] = []
    assume_fixed: dict[str, int] = {}
    for ts in succession:
        control_strategies.append(find_internal_drivers(bn, ts, assume_fixed))
        assume_fixed.update(ts)

    return control_strategies


def find_internal_drivers(
    bn: BooleanNetwork,
    target_trap_space: dict[str, int],
    assume_fixed: dict[str, int] | None = None,
) -> list[dict[str, int]]:
    """Finds internal drives of a given target trap space

    Parameters
    ----------
    bn : BooleanNetwork
        The network to analyze, which contains the Boolean update functions.
    target_trap_space : dict[str, int]
        The trap space we want to find drivers for.
    assume_fixed: dict[str,int] | None
        A dictionary of fixed variables that should be assumed to be fixed.


    Returns
    -------
    list[dict[str, int]]
        A list of internal driver sets, represented as state dictionaries.
    """
    if assume_fixed is None:
        assume_fixed = {}

    target_trap_space_inner = {
        k: v for k, v in target_trap_space.items() if k not in assume_fixed
    }
    max_drivers = len(target_trap_space_inner)
    drivers: list[dict[str, int]] = []
    for driver_set_size in range(1, max_drivers + 1):
        for driver_set in combinations(target_trap_space_inner, driver_set_size):
            if any(set(d) <= set(driver_set) for d in drivers):
                continue

            driver_dict = {k: target_trap_space_inner[k] for k in driver_set}
            ldoi, _ = percolate_space(bn, driver_dict | assume_fixed)

            if target_trap_space_inner.items() <= ldoi.items():
                drivers.append(driver_dict)

    return drivers
