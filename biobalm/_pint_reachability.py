"""
A module which provides utility methods for reachability analysis using Pint.

Note that in biobalm, we try to treat pint as an *optional* dependency. Hence, if you
want to use any of these methods, make sure to import the module in a way that fails
safely when Pint is not installed.
"""

from __future__ import annotations

from functools import reduce
from typing import TYPE_CHECKING, cast, Iterable

from networkx import DiGraph  # type: ignore
from pypint import Goal, InMemoryModel  # type:ignore

from biobalm.petri_net_translation import (
    place_to_variable,
    optimized_recursive_dnf_generator,
)
from biobalm.types import BooleanSpace, SuccessionDiagramConfiguration

if TYPE_CHECKING:
    from biodivine_aeon import Bdd


def pint_reachability(
    petri_net: DiGraph,
    initial_state: BooleanSpace,
    target_states: Bdd,
    config: SuccessionDiagramConfiguration,
) -> bool:
    """
    Use Pint to check if a given `initial_state` can possibly reach some state
    in the `target_states` BDD.

    If the reachability analysis is inconclusive, the method
    raises a `RuntimeError`.
    """
    if target_states.is_false():
        return False  # Cannot reach a stat in an empty set.

    # Build a Pint model through an automata network and copy
    # over the initial condition.
    pint_model = InMemoryModel(_petri_net_as_automata_network(petri_net))
    for var, level in initial_state.items():
        pint_model.initial_state[var] = level

    goal = _pint_build_symbolic_goal(target_states, config)

    def failed(*_):
        raise RuntimeError("Cannot verify.")

    return pint_model.reachability(goal=goal, fallback=failed)  # type: ignore


def _pint_build_symbolic_goal(
    states: Bdd, config: SuccessionDiagramConfiguration
) -> Goal:
    """
    A helper method which (very explicitly) converts a set of states
    represented through a BDD into a Pint `Goal`.

    Note that if `GOAL_SIZE_LIMIT` is exceeded, a partial goal is returned
    that may not cover all the states in the argument `Bdd`. This avoids
    exceeding the argument list size limit, but introduces additional source
    of incompleteness into the reachability process.
    """
    assert not states.is_false()

    goals: list[Goal] = []
    limit = config["pint_goal_size_limit"]
    for clause in optimized_recursive_dnf_generator(states):
        named_clause = {
            states.__ctx__().get_variable_name(var): int(value)
            for var, value in clause.items()
        }

        limit -= len(named_clause)
        if limit < 0:
            # If the goal is too large to be passed as a command line argument,
            # break here and don't continue. This is not ideal but I'm not sure
            # how to fix this other than modifying `pint` itself.
            if config["debug"]:
                print(
                    "WARNING: `pint` goal size limit exceeded. A partial goal is used."
                )
            break

        goal_atoms = [f"{var}={level}" for var, level in named_clause.items()]
        goals.append(Goal(",".join(goal_atoms)))

    return reduce(lambda a, b: a | b, goals)


def _petri_net_as_automata_network(petri_net: DiGraph) -> str:
    """
    Takes a Petri net which was created by implicant encoding from a Boolean network,
    and builds an automata network file (`.an`) compatible with the Pint tool.
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

        predecessors = set(cast(Iterable[str], petri_net.predecessors(transition)))  # type: ignore
        successors = set(cast(Iterable[str], petri_net.successors(transition)))  # type: ignore

        # The value under modification is the only
        # value that is different between successors and predecessors.
        source_place = next(iter(predecessors - successors))
        target_place = next(iter(successors - predecessors))

        (s_var, s_level) = place_to_variable(source_place)
        (t_var, t_level) = place_to_variable(target_place)
        assert s_var == t_var

        # The remaining places represent the necessary conditions.
        # Here, we transform them into a text format.
        condition_places = sorted(predecessors.intersection(successors))
        condition_tuples = [place_to_variable(p) for p in condition_places]
        conditions = [f'"{var}"={int(level)}' for var, level in condition_tuples]

        # A pint rule consists of a variable name, value transition,
        # and a list of necessary conditions for the transition (if any).
        if len(conditions) == 0:
            rule = f'"{s_var}" {int(s_level)} -> {int(t_level)}\n'
        else:
            rule = f"\"{s_var}\" {int(s_level)} -> {int(t_level)} when {' and '.join(conditions)}\n"

        automata_network += rule

    return automata_network
