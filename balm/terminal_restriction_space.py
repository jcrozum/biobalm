from __future__ import annotations

from typing import TYPE_CHECKING, cast

if TYPE_CHECKING:
    from biodivine_aeon import BooleanNetwork
    from pyeda.boolalg.bdd import BinaryDecisionDiagram

from pyeda.boolalg.bdd import expr2bdd
from pyeda.boolalg.expr import Not

from balm.drivers import find_single_drivers, find_single_node_LDOIs
from balm.pyeda_utils import aeon_to_pyeda
from balm.space_utils import percolate_network, percolation_conflicts
from balm.state_utils import state_list_to_bdd, state_to_bdd
from balm.trappist_core import trappist
from balm.types import space_type


def get_self_neg_tr_trap_spaces(network: BooleanNetwork) -> list[space_type]:
    """
    Takes a Boolean network and gets its self-negating time-reversal trap spaces.
    To find time-reversal trap spaces in a specific trap space,
    percolated network should be given as input.
    """
    tr_trap_spaces = trappist(network, problem="max", reverse_time=True)
    self_neg_tr_trap_spaces: list[space_type] = []
    for tr_trap_space in tr_trap_spaces:
        conflicts = percolation_conflicts(network, tr_trap_space)
        if conflicts:
            self_neg_tr_trap_spaces.append(tr_trap_space)

    return self_neg_tr_trap_spaces


def get_terminal_restriction_space(
    stable_motifs: list[space_type],
    network: BooleanNetwork,
    ensure_subspace: space_type,
    use_single_node_drivers: bool = True,
    use_tr_trapspaces: bool = True,
) -> BinaryDecisionDiagram:
    """
    Find the terminal restriction space.

    Construct the negation of the terminal restriction space and return the negation of the negation.
    ~terminal restriction space = stable motifs | ~R(X) | self negating time reversal trap spaces
    ~R(X) = single node drivers of stable motifs(delta) | F(delta) | ~LDOI(~delta)

    network and ensure subspace required for using delta and tr trapspaces.
    """

    # ~terminal restriction space includes stable motifs
    result_bdd = state_list_to_bdd(stable_motifs)

    if use_single_node_drivers or use_tr_trapspaces:
        # Get the percolated Boolean Network
        reduced_network = percolate_network(network, ensure_subspace)

        # ~terminal restriction space includes ~R(X)
        if use_single_node_drivers:
            LDOIs = find_single_node_LDOIs(reduced_network)

            for stable_motif in stable_motifs:
                # find delta
                single_node_drivers = find_single_drivers(
                    stable_motif, reduced_network, LDOIs=LDOIs
                )
                if len(single_node_drivers) == 0:
                    continue

                # ~R(X) includes delta
                result_bdd = result_bdd | state_to_bdd(
                    cast(space_type, dict(single_node_drivers))
                )

                for single_node_driver in single_node_drivers:
                    # ~R(X) includes the F(delta)
                    expression = aeon_to_pyeda(
                        reduced_network.get_update_function(single_node_driver[0])
                    )
                    if single_node_driver[1] == 0:
                        expression = Not(expression)
                    result_bdd = result_bdd | expr2bdd(expression)

                    # get ~delta
                    not_delta: tuple[str, int] = (
                        single_node_driver[0],
                        1 - single_node_driver[1],
                    )
                    # ~R(X) includes ~LDOI(~delta)
                    result_bdd = result_bdd | ~state_to_bdd(LDOIs[not_delta])

        # ~terminal restriction space includes self negating time reversal trapspaces
        if use_tr_trapspaces:
            self_neg_tr_trap_spaces = get_self_neg_tr_trap_spaces(reduced_network)
            if len(self_neg_tr_trap_spaces) == 0:
                return ~result_bdd
            else:
                result_bdd = result_bdd | state_list_to_bdd(self_neg_tr_trap_spaces)

    return ~result_bdd