import random # importing the random module

from pyeda.boolalg import boolfunc
from pyeda.boolalg.bdd import bddvar, expr2bdd
from pyeda.boolalg.expr import expr

from typing import List, Set

def motif_avoidant_check(candidates, all_traps, U_neg: List[str]):
    B = set_retained_set(U_neg, candidates)

    # for node in B.keys():
    #     print(node + " = " + str(B[node]))

    F = compute_candidate_set(U_neg, B, candidates)

    if not F.is_empty():
        F = PreprocessingSSF(F, all_traps)

        if not F.is_empty():
            # print ("Start the filtering process")

            A = FilteringProcess(F, all_traps)


def set_retained_set(U_neg: List[str], candidates):
    B = {}

    for node in U_neg:
        B[node] = random.randint(0, 1)

        """
        TODO: More heuristics for setting the retained set
        """

    return B


def compute_candidate_set(U_neg: List[str], B, candidates):
    F = candidates

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
