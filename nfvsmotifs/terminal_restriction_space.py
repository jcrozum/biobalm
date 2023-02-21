from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from biodivine_aeon import BooleanNetwork # type: ignore    

from nfvsmotifs.trappist_core import trappist
from nfvsmotifs.space_utils import percolate_space


def get_self_neg_tr_trap_spaces(network: BooleanNetwork) -> list[dict[str, int]]:
    """
    Takes a Boolean network and gets its self-negating time-reversal trap spaces.
    To find time-reversal trap spaces in a specific trap space,
    percolated network should be given as input. 
    """
    tr_trap_spaces = trappist(network, problem="max", reverse_time=True)
    self_neg_tr_trap_spaces = []
    for tr_trap_space in tr_trap_spaces:
        result, conflicts = percolate_space(network, tr_trap_space)
        if conflicts:
            self_neg_tr_trap_spaces.append(tr_trap_space)

    return self_neg_tr_trap_spaces
