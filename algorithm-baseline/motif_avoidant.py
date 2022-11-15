import random # importing the random module

from pyeda.boolalg import boolfunc
from pyeda.boolalg.bdd import bddvar, expr2bdd
from pyeda.boolalg.expr import expr

from typing import List, Set

from trappist import compute_trap_spaces

def motif_avoidant_check(candidates, all_traps, U_neg: List[str], petri_net, subspace, source_nodes):
    trs = get_terminal_restriction_space(petri_net, subspace, all_traps, source_nodes)

    B = set_retained_set(U_neg, trs)

    # for node in B.keys():
    #     print(node + " = " + str(B[node]))

    F = compute_candidate_set(U_neg, B, trs)

    if not F.is_empty():
        F = PreprocessingSSF(F, all_traps)

        if not F.is_empty():
            # print ("Start the filtering process")

            A = FilteringProcess(F, all_traps)

def get_terminal_restriction_space(petri_net, subspace, all_traps, source_nodes):
    trs = all_traps

    tr_trap_spaces = get_tr_trap_spaces(petri_net, subspace, source_nodes)
    self_neg_tr_trap_spaces = get_self_neg_tr_trap_spaces(tr_trap_spaces)

    delta = find_single_drivers(all_traps)

    return trs


def get_tr_trap_spaces(petri_net, subspace, source_nodes):
    """
    TODO: get time-reversal trap spaces
    """

    tr_trap_spaces = []

    #tr_trap_spaces = compute_trap_spaces(petri_net, computation="max", subspace=subspace, source_nodes=source_nodes, time_reversal=True)

    # for trap in tr_trap_spaces:
    #     print(trap)

    return tr_trap_spaces


def get_self_neg_tr_trap_spaces(tr_trap_spaces):
    """
    TODO: get self-negating time-reversal trap spaces by percolating
    """

    self_neg_tr_trap_spaces = []

    return self_neg_tr_trap_spaces

def find_single_drivers(all_traps):
    """
    TODO: find single drivers of trap spaces
    """

    delta = []

    return delta


def set_retained_set(U_neg: List[str], trs):
    B = {}

    for node in U_neg:
        B[node] = random.randint(0, 1)

        """
        TODO: More heuristics for setting the retained set
        """

    return B


def compute_candidate_set(U_neg: List[str], B, trs):
    F = trs

    """
    TODO: compute the set of fixed points of the reduced STG lying on candidates
    """

    return F


def PreprocessingSSF(F, all_traps):
    F_new = F

    """
    TODO: Preprocessing SSF with the target set all_traps
    """

    return F_new


def FilteringProcess(F, all_traps):
    A = [] # the set of motif-avoidant attractors

    """
    TODO: Filtering out the candidate set by using reachability analysis
    """

    return A
