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

    if len(F) > 0:
        F = PreprocessingSSF(network, F, terminal_res_space)

        if len(F) > 0:
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


def FilteringProcess(network: BooleanNetwork, petri_net: DiGraph, F: list[dict[str, int]], terminal_res_space: BinaryDecisionDiagram) -> list[dict[str, int]]:
    list_motif_avoidant_atts: list[dict[str, int]] = []

    """
        TODO: Filtering out the candidate set by using the reachability analysis.
    """

    target_set = ~terminal_res_space
    F_bdd = list_state_2_bdd(F)
    A = 0

    while len(F) > 0:
        state = F.pop(0)

        state_bdd = state_2_bdd(state)
        F_bdd = F_bdd & ~state_bdd

        joint_target_set = target_set | F_bdd | A

        if ABNReach(network, petri_net, state, joint_target_set) == False:
            A = A | state_bdd
            list_motif_avoidant_atts.append(state)
        

    return list_motif_avoidant_atts


def ABNReach(network: BooleanNetwork, petri_net: DiGraph, state: dict[str, int], joint_target_set: BinaryDecisionDiagram) -> bool:
    is_reachable: bool = False

    # The first phase using Pint <https://loicpauleve.name/pint/doc/transient-analysis.html>
    pint_result = PintReach(network, state, joint_target_set)

    if pint_result == "True":
        is_reachable = True
    elif pint_result == "False":
        is_reachable = False
    else:
        # The second phase using SAT-based bounded model checking
        d_bound: int = 20
        sat_result = SATReach(network, state, joint_target_set, d_bound)

        if sat_result == "True":
            is_reachable = True
        else:
            # The final phase using Petri net unfoldings
            # This phase is the last resort ensuring the correctness of ABNReach

            is_reachable = MoleReach(network, petri_net, state, joint_target_set)

    return is_reachable


def PintReach(network: BooleanNetwork, state: dict[str, int], joint_target_set: BinaryDecisionDiagram) -> str:
    # TODO

    return "False"


def SATReach(network: BooleanNetwork, state: dict[str, int], joint_target_set: BinaryDecisionDiagram, d_bound: int) -> str:
    # TODO

    return "False"


def MoleReach(network: BooleanNetwork, petri_net: DiGraph, state: dict[str, int], joint_target_set: BinaryDecisionDiagram) -> bool:
    # TODO
    
    return False
