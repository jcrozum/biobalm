import random # type: ignore

from pyeda.boolalg import boolfunc # type:ignore
from pyeda.boolalg.bdd import bddvar, expr2bdd, BinaryDecisionDiagram # type:ignore
from pyeda.boolalg.expr import expr # type:ignore

from biodivine_aeon import BooleanNetwork # type: ignore

from typing import List, Set, Dict # type: ignore

from networkx import DiGraph # type: ignore

from nfvsmotifs.pyeda_utils import aeon_to_pyeda # type:ignore
from nfvsmotifs.state_utils import state_2_bdd, list_state_2_bdd, eval_function, is_member_bdd # type:ignore

"""
    A state is represented as a dict.
    The considered Boolean network is represented by a BooleanNetwork object provided by AEON.
    The terminal restriction space is represented as a BDD.
"""

def motif_avoidant_check(network: BooleanNetwork, petri_net: DiGraph, F: list[dict[str, int]], terminal_res_space: BinaryDecisionDiagram) -> list[dict[str, int]]:
    """
        Return the list of states corresponding to motif-avoidant attractors.
        This list may be empty, indicating that there are no motif-avoidant attractors.
    """
    list_motif_avoidant_atts = []

    if not F.is_empty():
        F = PreprocessingSSF(network, F, terminal_res_space)

        if not F.is_empty():
            """
                PreprocessingSSF does not reach the best case.
                Hence we need to use the reachability analysis on the asynchronous Boolean network.
            """

            list_motif_avoidant_atts = FilteringProcess(network, petri_net, F, terminal_res_space)


    return list_motif_avoidant_atts


def PreprocessingSSF(network: BooleanNetwork, F: list[dict[str, int]], terminal_res_space: BinaryDecisionDiagram) -> list[dict[str, int]]:
    F_result = []
    I_MAX = 2

    target_set = ~terminal_res_space
    F_bdd = list_state_2_bdd(F)

    nodes = []
    funs_bdd = {}
    for var in network.variables():
        var_name = network.get_variable_name(var)
        nodes.append(var_name)

        function = network.get_update_function(var)
        function = aeon_to_pyeda(function)

        fx = expr2bdd(function)
        funs_bdd[var_name] = fx


    for state in F:
        state_bdd = state_2_bdd(state)
        F_bdd = F_bdd & ~state_bdd

        reach_target_set = False
        for i in range(1, I_MAX + 1):
            random.shuffle(nodes)

            for node in nodes:
                #print(node)
                state[node] = eval_function(funs_bdd[node], state)

            #print(state)

            if is_member_bdd(state, F_bdd) or is_member_bdd(state, target_set):
                reach_target_set = True
                #print("Reach terget set")
                break

        if reach_target_set == False:
            F_bdd = F_bdd | state_bdd
            F_result.append(state)


    return F_result


def FilteringProcess(network: BooleanNetwork, petri_net: DiGraph, F: list[dict[str, int]], terminal_res_space: BinaryDecisionDiagram):
    list_motif_avoidant_atts = []

    """
        TODO: Filtering out the candidate set by using the reachability analysis.
    """

    return list_motif_avoidant_atts
