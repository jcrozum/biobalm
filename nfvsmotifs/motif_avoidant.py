import random # type: ignore
import os # type: ignore

from pyeda.boolalg import boolfunc # type:ignore
from pyeda.boolalg.bdd import bddvar, expr2bdd, BinaryDecisionDiagram # type:ignore
from pyeda.boolalg.expr import expr # type:ignore

from biodivine_aeon import BooleanNetwork # type: ignore

from typing import List, Set, Dict # type: ignore

from networkx import DiGraph # type: ignore

from nfvsmotifs.pyeda_utils import aeon_to_pyeda # type:ignore
from nfvsmotifs.state_utils import state_2_bdd, list_state_2_bdd, eval_function, is_member_bdd # type:ignore

import pypint # type:ignore
from pypint import Goal # type:ignore

"""
    A state is represented as a dict.
    The considered Boolean network is represented by a BooleanNetwork object provided by AEON.
    The terminal restriction space is represented as a BDD.
"""

def motif_avoidant_check(network: BooleanNetwork, petri_net: DiGraph, F: list[dict[str, int]], terminal_res_space: BinaryDecisionDiagram, bn_name: str) -> list[dict[str, int]]:
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

            list_motif_avoidant_atts = FilteringProcess(network, petri_net, F, terminal_res_space, bn_name)


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
                state[node] = eval_function(funs_bdd[node], state)

            if is_member_bdd(state, F_bdd) or is_member_bdd(state, target_set):
                reach_target_set = True
                break

        if reach_target_set == False:
            F_bdd = F_bdd | state_bdd
            F_result.append(state)


    return F_result


def FilteringProcess(network: BooleanNetwork, petri_net: DiGraph, F: list[dict[str, int]], terminal_res_space: BinaryDecisionDiagram, bn_name: str) -> list[dict[str, int]]:
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

        if ABNReach_current_version(network, state, joint_target_set, bn_name) == False:
            A = A | state_bdd
            list_motif_avoidant_atts.append(state)
        

    return list_motif_avoidant_atts


def Pint_get_goal_from_bdd(joint_target_set: BinaryDecisionDiagram) -> Goal:
    list_path = []

    for i, t in enumerate(joint_target_set.satisfy_all()):
        list_path.append(t)


    if len(list_path) > 0:
        t = list_path.pop(0)

        node_goal_str = []

        for p, v in t.items():
            if v == 0:
                node_goal_str.append(str(p) + "=0")
            else:
                node_goal_str.append(str(p) + "=1")

        goal = Goal(",".join(node_goal_str))


    for t in list_path:
        node_goal_str = []

        for p, v in t.items():
            if v == 0:
                node_goal_str.append(str(p) + "=0")
            else:
                node_goal_str.append(str(p) + "=1")

        goal = goal | Goal(",".join(node_goal_str))


    return goal


def ABNReach_current_version(network: BooleanNetwork, state: dict[str, int], joint_target_set: BinaryDecisionDiagram, bn_name: str) -> bool:
    # In this version, we only use Pint.
    # If the result of Pint is True (i.e., reachable) or False (i.e., unreachable), we simply return this result.
    # If the result of Pint is Inconc (i.e., not determining whether reachable or unreachable), Pint falls back to an exact model checker (e.g., Mole).
    # However, in the latter case, we can exploit the goal reduction technique of Pint.

    is_reachable: bool = False

    # Assume that we have the .an file in the current dictory
    # In the future, we can build the asynchronous automaton (input of Pint) from the BooleanNetwork object)
    an_file = os.getcwd() + "/" + bn_name + ".an"

    m = pypint.load(an_file)

    # set the initial state
    for node in state:
        m.initial_state[node] = state[node]

    # get the goal
    goal = Pint_get_goal_from_bdd(joint_target_set)

    is_reachable = m.reachability(goal=goal, fallback='mole')

    return is_reachable


def ABNReach_under_review_version(network: BooleanNetwork, petri_net: DiGraph, state: dict[str, int], joint_target_set: BinaryDecisionDiagram, bn_name: str) -> bool:
    is_reachable: bool = False

    # The first phase using Pint <https://loicpauleve.name/pint/doc/transient-analysis.html>
    pint_result = PintReach(network, state, joint_target_set, bn_name)

    if pint_result == "True":
        is_reachable = True
    elif pint_result == "False":
        is_reachable = False
    else:
        # The second phase using SAT-based bounded model checking
        d_bound: int = 20
        sat_result = SATReach(network, state, joint_target_set, d_bound, bn_name)

        if sat_result == "True":
            is_reachable = True
        else:
            # The final phase using Petri net unfoldings
            # This phase is the last resort ensuring the correctness of ABNReach

            is_reachable = MoleReach(network, petri_net, state, joint_target_set, bn_name)

    return is_reachable


def PintReach(network: BooleanNetwork, state: dict[str, int], joint_target_set: BinaryDecisionDiagram, bn_name: str) -> str:
    # Install Pint on Linux by: `conda install -c colomoto pint`
    # Install the Python interface of Pint by: `pip3 install -U pypint`

    # Assume that we have the .an file in the current dictory
    # In the future, we can build the asynchronous automaton (input of Pint) from the BooleanNetwork object)
    an_file = os.getcwd() + "/" + bn_name + ".an"

    m = pypint.load(an_file)

    # set the initial state
    for node in state:
        m.initial_state[node] = state[node]

    # get the goal
    goal = Pint_get_goal_from_bdd(joint_target_set)

    pint_result = m.reachability(goal=goal, fallback=None)

    if pint_result == pypint.types.Inconc:
        pint_result = "Inconc"
    elif pint_result == True:
        pint_result = "True"
    else:
        pint_result = "False"

    return pint_result


def SATReach(network: BooleanNetwork, state: dict[str, int], joint_target_set: BinaryDecisionDiagram, d_bound: int, bn_name: str) -> str:
    # TODO

    return "False"


def MoleReach(network: BooleanNetwork, petri_net: DiGraph, state: dict[str, int], joint_target_set: BinaryDecisionDiagram, bn_name: str) -> bool:
    # TODO
    
    return False
