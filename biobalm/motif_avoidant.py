"""
A module which provides utility methods for detecting motif-avoidant attractors
within a terminal restriction space.
"""

from __future__ import annotations

import random
from functools import reduce
from typing import TYPE_CHECKING

from networkx import DiGraph  # type: ignore
from pypint import Goal, InMemoryModel  # type:ignore

from biobalm.petri_net_translation import place_to_variable
from biobalm.space_utils import dnf_function_is_true, remove_state_from_dnf
from biobalm.symbolic_utils import (
    function_eval,
    function_is_true,
    state_list_to_bdd,
    state_to_bdd,
)
from biobalm.types import BooleanSpace

if TYPE_CHECKING:
    from biodivine_aeon import AsynchronousGraph, Bdd


def make_retained_set(
    graph: AsynchronousGraph,
    nfvs: list[str],
    space: BooleanSpace,
    child_spaces: list[BooleanSpace] | None = None,
) -> BooleanSpace:
    """
    Calculate the retained set for a given `space`.

    The retained set is technically a space-like object that describes the
    variables which have to be fixed in order for the network to lose all
    complex attractors. However, note that this really means changing the update
    functions. This is not a trap space that only contains fixed-points,
    but a description of how the network must be modified to eliminate
    all complex attractors in the given `space`.

    The construction guarantees that any complex attractor of the old
    network will manifest as at least one fixed-point in the new network. But this
    network modification is not implemented here. This method only generates
    the necessary list of variables and values.

    Finally, note that the method uses a heuristic to select values
    that should lead to the least amount of fixed-points in
    the modified network.

    Parameters
    ----------
    graph : AsynchronousGraph
        The symbolic update functions stored as an `AsynchronousGraph` object
        from the `biodivine_aeon` library.
    nfvs : list[str]
        The list of variables in the NFVS that is valid for a network
        restricted to the given `space`.
    space : BooleanSpace
        A :class:`BooleanSpace<biobalm.types.BooleanSpace>` object describing the
        current trap space.
    child_spaces : list[BooleanSpace] | None, optional
        A list of :class:`BooleanSpace<biobalm.types.BooleanSpace>` objects
        describing the child spaces of the current node. Only attractors that
        are not in the child spaces are considered. If no child spaces are
        provided, then all attractors are considered.

    Returns
    -------
    BooleanSpace
        A :class:`BooleanSpace<biobalm.types.BooleanSpace>` object describing the
        retained set.
    """

    if child_spaces is None:
        child_spaces = []

    # Initially, the retained set only contains the fixed values from the
    # current node space (this elimiantes unnecessary Petri net transitions for
    # values which we already proved are constant).
    #
    # In the following code, we then extend the retained set based on the
    # model's NFVS and the current child spaces.
    retained_set = space.copy()

    # First, if there are any child spaces present, we extend the retained set
    # with the values from the one that has the least amount of fixed variables
    # shared with the NFVS.
    if len(child_spaces) > 0:
        # Find the child space that has the fewest nodes in common with the NFVS:
        least_common_child_space = child_spaces[0]
        least_common_nodes = len(set(least_common_child_space) & set(nfvs))
        for child_space in child_spaces:
            common_nodes = len(set(child_space) & set(nfvs))
            if common_nodes < least_common_nodes:
                least_common_nodes = common_nodes
                least_common_child_space = child_space

        for x in least_common_child_space:
            if (x not in retained_set) and (x in nfvs):
                retained_set[x] = least_common_child_space[x]

    # Then, set the remaining NFVS variables based on the majority output value
    # in the update function of the relevant variable.
    for x in nfvs:
        if x in retained_set:
            continue

        fn_bdd = graph.mk_update_function(x)

        # If most of the function is positive, we set the value
        # to `1`, otherwise set it to `0`.
        if fn_bdd.cardinality() > fn_bdd.l_not().cardinality():
            retained_set[x] = 1
        else:
            retained_set[x] = 0

    return retained_set


def detect_motif_avoidant_attractors(
    graph: AsynchronousGraph,
    petri_net: DiGraph,
    candidates: list[BooleanSpace],
    terminal_restriction_space: Bdd,
    max_iterations: int,
    is_in_an_mts: bool = False,
) -> list[BooleanSpace]:
    """
    Determine which seed states from the `candidates` list map to true
    motif-avoidant attractors.

    Assumes that all attractors in the `terminal_restriction_space` are
    covered by at least one state in `candidates`. Eliminates all states
    from `candidates` that are provably not in an attractor, or that correspond
    to an attractor that is already covered by some other candidate.

    The method consists of two steps: initial pruning, and exact reachability
    analysis. As such, the result is always complete. However, some of the
    parameters can be used to adjust the initial pruning step.

    If `is_in_an_mts` is set, the method assumes the `candidates` are members
    of a minimal trap space. In such case, a different simulation algorithm
    is used for initial pruning. Furthermore, `max_iterations` can be used to
    to limit the time spent on the initial pruning step.

    Parameters
    ----------
    graph : AsynchronousGraph
        The symbolic update functions stored as an `AsynchronousGraph` object
        from the `biodivine_aeon` library.
    petri_net : DiGraph
        The Petri net representation of the update functions.
    candidates : list[BooleanSpace]
        A list of :class:`BooleanSpace<biobalm.types.BooleanSpace>` objects
        describing the candidate seed states.
    terminal_restriction_space : Bdd
        A symbolic set of states which contains all motif avoidant attractors
        in question (i.e. if a candidate state can leave this set, the candidate
        cannot be an attractor).
    max_iterations : int
        Specifies how much time should be spent on the simpler preprocessing methods.
    is_in_an_mts : bool, optional
        By default `False`. If `True`, then we assume that the candidates lie
        within a minimal trap space, enabling certain optimizations. If
        incorrectly set to `True`, the result is undefined.

    Returns
    -------
    list[BooleanSpace]
        A list of :class:`BooleanSpace<biobalm.types.BooleanSpace>` objects
        describing the motif-avoidant attractors among the input candidate set.
    """
    if len(candidates) == 0:
        return []

    if len(candidates) == 1 and is_in_an_mts:
        return candidates

    candidates = _preprocess_candidates(
        graph,
        candidates,
        terminal_restriction_space,
        max_iterations,
        is_in_an_mts=is_in_an_mts,
    )

    if len(candidates) == 0:
        return []

    if len(candidates) == 1 and is_in_an_mts:
        return candidates

    return _filter_candidates(petri_net, candidates, terminal_restriction_space)


def _preprocess_candidates(
    graph: AsynchronousGraph,
    candidates: list[BooleanSpace],
    terminal_restriction_space: Bdd,
    max_iterations: int,
    is_in_an_mts: bool = False,
    simulation_seed: int = 0,
) -> list[BooleanSpace]:
    """
    A fast but incomplete method for eliminating spurious attractor candidates.

    The idea is to build the symbolic encoding of the given `network`, and then
    randomly simulate transitions for individual states, trying to reach a state
    outside of the `terminal_restriction_space`.

    TODO (2): We could probably make this algorithm slighlty less random by doing
    a limited version of symbolic reachability. I.e. instead of simulating just one
    state transition in each step, compute the whole successor BDD and then test
    against that. Once the BDD becomes too large after several steps, we can just
    pick a single state from it and start again. Sam: I'll add a version of this
    later, once we can actually benchmark how it performs :)
    """
    # A random generator initialized with a fixed seed. Ensures simulation
    # is randomized but deterministic.
    generator = random.Random(simulation_seed)

    variables = graph.network_variable_names()
    update_functions = {var: graph.mk_update_function(var) for var in variables}

    # Use stochastic simulation to prune the set of candidate states.
    # We use different simulation approach depending on whether this space
    # is a minimal trap or not. In previous work, this was shown to work
    # well, but in the future we need to better document the resoning
    # behind these two algorithms.
    if not is_in_an_mts:
        # Copy is sufficient because we won't be modifying the states within the set.
        candidates_dnf = candidates.copy()
        filtered_candidates: list[BooleanSpace] = []
        for state in candidates:
            # Remove the state from the candidates. If we can prove that is
            # is not an attractor, we will put it back.
            candidates_dnf = remove_state_from_dnf(candidates_dnf, state)

            simulation = state.copy()  # A copy of the state that we can overwrite.
            is_valid_candidate = True
            for _ in range(max_iterations):
                # Advance all variables by one step in random order.
                generator.shuffle(variables)
                for var in variables:
                    step = function_eval(update_functions[var], simulation)
                    assert step is not None
                    simulation[var] = step

                if dnf_function_is_true(candidates_dnf, simulation):
                    # The state can reach some other state in the candidate
                    # set. This does not mean it cannot be an attractor, but
                    # it means it is sufficient to keep considering
                    # the remaining candidates.
                    is_valid_candidate = False
                    break

                if not function_is_true(terminal_restriction_space, simulation):
                    # The state can reach some other state outside of the
                    # terminal restriction space, which means it cannot be
                    # a motif avoidant attractor in this subspace.
                    is_valid_candidate = False
                    break

            if is_valid_candidate:
                # If we cannot rule out the candidate, we have to put it back
                # into the candidate set.
                candidates_dnf.append(state)
                filtered_candidates.append(state)

        return filtered_candidates
    else:
        filtered_candidates = []
        for _ in range(max_iterations):
            generator.shuffle(variables)
            candidates_dnf = candidates.copy()
            filtered_candidates = []

            for state in candidates:
                candidates_dnf = remove_state_from_dnf(candidates_dnf, state)

                simulation = state.copy()
                for var in variables:
                    step = function_eval(update_functions[var], simulation)
                    assert step is not None
                    simulation[var] = step

                if not dnf_function_is_true(candidates_dnf, simulation):
                    candidates_dnf.append(simulation)
                    filtered_candidates.append(simulation)

            if len(filtered_candidates) <= 1:
                break

            candidates = filtered_candidates

        return filtered_candidates


def _filter_candidates(
    petri_net: DiGraph,
    candidates: list[BooleanSpace],
    terminal_restriction_space: Bdd,
) -> list[BooleanSpace]:
    """
    Filter candidate states using reachability procedure in Pint.
    """

    ctx = terminal_restriction_space.__ctx__()
    candidates_bdd = state_list_to_bdd(ctx, candidates)
    avoid_states = candidates_bdd.l_and_not(terminal_restriction_space)
    filtered_candidates: list[BooleanSpace] = []

    for state in candidates:
        state_bdd = state_to_bdd(ctx, state)

        # Remove state from the symbolic set. If we can prove that it
        # is not an attractor, we will put it back.
        avoid_states = avoid_states.l_and_not(state_bdd)

        if not _Pint_reachability(petri_net, state, avoid_states):
            avoid_states = avoid_states.l_or(state_bdd)
            filtered_candidates.append(state)

    return filtered_candidates


def _Pint_reachability(
    petri_net: DiGraph,
    initial_state: BooleanSpace,
    target_states: Bdd,
) -> bool:
    """
    Use Pint to check if a given `initial_state` can possibly reach some state
    in the `target_states` BDD.

    TODO: Here, if the result of static analysis is inconclusive, Pint falls back to `mole`
    model checker. However, in the future, we might also explore other methods, such as
    petri net reduction or symbolic reachability.
    """
    if target_states.is_false():
        return False  # Cannot reach a stat in an empty set.

    # Build a Pint model through an automata network and copy
    # over the initial condition.
    pint_model = InMemoryModel(_petri_net_as_automata_network(petri_net))
    for var, level in initial_state.items():
        pint_model.initial_state[var] = level

    goal = _Pint_build_symbolic_goal(target_states)

    return pint_model.reachability(goal=goal, fallback="mole")  # type: ignore


def _Pint_build_symbolic_goal(states: Bdd) -> Goal:
    """
    A helper method which (very explicitly) converts a set of states
    represented through a BDD into a Pint `Goal`.
    """
    assert not states.is_false()

    goals: list[Goal] = []
    for clause in states.clause_iterator():
        named_clause = {
            states.__ctx__().get_variable_name(var): int(value)
            for var, value in clause.items()
        }
        goal_atoms = [f"{var}={level}" for var, level in named_clause.items()]
        goals.append(Goal(",".join(goal_atoms)))

    return reduce(lambda a, b: a | b, goals)


def _petri_net_as_automata_network(petri_net: DiGraph) -> str:
    """
    Takes a Petri net which was created by implicant encoding from a Boolean network,
    and builds an automata network file (`.an`) compatible with the Pint tool.

    TODO: This is one of those things that would probably be better served by having
    an "explicit" `PetriNetEncoding` class.
    """
    automata_network = ""

    # Go through all PN places and save them as model variables.
    variable_set: set[str] = set()
    for place, kind in petri_net.nodes(data="kind"):  # type: ignore
        if kind != "place":
            continue
        variable_set.add(place_to_variable(place)[0])  # type: ignore[reportUnknownArgumentType] # noqa
    variables = sorted(variable_set)

    # Declare all variables with 0/1 domains.
    for var in variables:
        automata_network += f'"{var}" [0, 1]\n'

    for transition, kind in petri_net.nodes(data="kind"):  # type: ignore
        if kind != "transition":
            continue

        predecessors = set(petri_net.predecessors(transition))  # type: ignore
        successors = set(petri_net.successors(transition))  # type: ignore

        # The value under modification is the only
        # value that is different between successors and predecessors.
        source_place = next(iter(predecessors - successors))  # type: ignore[reportUnknownVariableType,reportUnknownArgumentType] # noqa
        target_place = next(iter(successors - predecessors))  # type: ignore[reportUnknownVariableType,reportUnknownArgumentType] # noqa

        (s_var, s_level) = place_to_variable(source_place)  # type: ignore[reportUnknownArgumentType] # noqa
        (t_var, t_level) = place_to_variable(target_place)  # type: ignore[reportUnknownArgumentType] # noqa
        assert s_var == t_var

        # The remaining places represent the necessary conditions.
        # Here, we transform them into a text format.
        conditions = sorted(predecessors.intersection(successors))  # type: ignore[reportUnknownVariableType,reportUnknownArgumentType] # noqa
        conditions = [place_to_variable(p) for p in conditions]  # type: ignore[reportUnknownVariableType,reportUnknownArgumentType] # noqa
        conditions = [f'"{var}"={int(level)}' for var, level in conditions]

        # A pint rule consists of a variable name, value transition,
        # and a list of necessary conditions for the transition (if any).
        if len(conditions) == 0:
            rule = f'"{s_var}" {int(s_level)} -> {int(t_level)}\n'
        else:
            rule = f"\"{s_var}\" {int(s_level)} -> {int(t_level)} when {' and '.join(conditions)}\n"

        automata_network += rule

    return automata_network
