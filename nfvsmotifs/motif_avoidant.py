import random # type: ignore
import os # type: ignore
import tempfile # type: ignore

from pyeda.boolalg import boolfunc # type:ignore
from pyeda.boolalg.bdd import bddvar, expr2bdd, BinaryDecisionDiagram # type:ignore
from pyeda.boolalg.expr import expr # type:ignore

from biodivine_aeon import BooleanNetwork # type: ignore

from typing import List, Set, Dict, IO # type: ignore

from networkx import DiGraph # type: ignore

from nfvsmotifs.pyeda_utils import aeon_to_pyeda # type:ignore
from nfvsmotifs.state_utils import state_2_bdd, list_state_2_bdd, eval_function, is_member_bdd # type:ignore
from nfvsmotifs.petri_net_translation import network_to_petrinet # type:ignore

import pypint # type:ignore
from pypint import Goal # type:ignore

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


def PreprocessingSSF(network: BooleanNetwork, F: list[dict[str, int]], terminal_res_space: BinaryDecisionDiagram, I_MAX: int) -> list[dict[str, int]]:
    """
        I_MAX: the maximum number of iterations of PreprocessingSSF
    """

    F_result = []

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


def FilteringProcess(network: BooleanNetwork, petri_net: DiGraph, F: list[dict[str, int]], terminal_res_space: BinaryDecisionDiagram) -> list[dict[str, int]]:
    list_motif_avoidant_atts: list[dict[str, int]] = []

    """
        Filtering out the candidate set by using the reachability analysis.
    """

    target_set = ~terminal_res_space
    F_bdd = list_state_2_bdd(F)
    A = 0

    while len(F) > 0:
        state = F.pop(0)

        state_bdd = state_2_bdd(state)
        F_bdd = F_bdd & ~state_bdd

        joint_target_set = target_set | F_bdd | A

        if ABNReach_current_version(network, petri_net, state, joint_target_set) == False:
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


def write_an_file_from_petri_net(network: BooleanNetwork, petri_net: DiGraph, an_file: IO):
    nodes = []
    for var in network.variables():
        var_name = network.get_variable_name(var)
        nodes.append(var_name)

    # write list of automata
    domain = "[0, 1]"
    for node in nodes:
        print(
            f"\"{node}\" {domain}", file=an_file
        )

    # write list of transitions
    for node, kind in petri_net.nodes(data="kind"):
        if kind == "place":
            continue
        else:  # it's a transition
            preds = list(petri_net.predecessors(node))
            succs = list(petri_net.successors(node))

            # value change
            source_places = list(set(preds) - set(succs))
            source_place = source_places[0]
            target_places = list(set(succs) - set(preds))
            target_place = target_places[0]

            # the list of conditions for the transition
            preds.remove(source_place)
            succs.remove(target_place)

            conds = []
            for pred in preds:
                node_name = "\"" + pred[3:] + "\""
                value = pred[1]
                cond = node_name + "=" + value

                conds.append(cond)

            conds = " and ".join(conds)

            value_change = source_place[1] + " -> " + target_place[1]

            # node name
            node_name = source_place[3:]

            # write the transition
            print(
                f"\"{node_name}\" {value_change} when {conds}", file=an_file
            )


def ABNReach_current_version(network: BooleanNetwork, petri_net: DiGraph, state: dict[str, int], joint_target_set: BinaryDecisionDiagram) -> bool:
    # In this version, we only use Pint.
    # If the result of Pint is True (i.e., reachable) or False (i.e., unreachable), we simply return this result.
    # If the result of Pint is Inconc (i.e., not determining whether reachable or unreachable), Pint falls back to an exact model checker (e.g., Mole).
    # However, in the latter case, we can exploit the goal reduction technique of Pint.

    is_reachable: bool = False
 
    (fd, tmpname) = tempfile.mkstemp(suffix=".an", text=True)
    with open(tmpname, "wt") as an_file:
        write_an_file_from_petri_net(network, petri_net, an_file)

    m = pypint.load(tmpname)

    # set the initial state
    for node in state:
        m.initial_state[node] = state[node]

    # get the goal
    goal = Pint_get_goal_from_bdd(joint_target_set)

    is_reachable = m.reachability(goal=goal, fallback='mole')

    # close .an file after finishing reachability analysis
    os.close(fd)
    os.unlink(tmpname)

    return is_reachable

