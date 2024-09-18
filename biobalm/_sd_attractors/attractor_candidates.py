"""
Provides functions that are used to identify and reduce the set of attractor candidate states.

Note that most of these functions assume they are operating on the network or petri net
percolated to the subspace of the relevant node.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal, cast

if TYPE_CHECKING:
    from biobalm.succession_diagram import SuccessionDiagram
    from biobalm.types import BooleanSpace
    from networkx import DiGraph  # type: ignore

import random
import biobalm
from biodivine_aeon import Bdd, AsynchronousGraph, BddVariable
from biobalm.trappist_core import compute_fixed_point_reduced_STG
from biobalm.symbolic_utils import state_list_to_bdd, valuation_to_state, state_to_bdd
from biobalm.space_utils import intersect, is_subspace

try:
    pint_available = True
    import biobalm._pint_reachability
except ModuleNotFoundError:
    pint_available = False


def compute_attractor_candidates(
    sd: SuccessionDiagram,
    node_id: int,
    greedy_asp_minification: bool,
    simulation_minification: bool,
    pint_minification: bool,
) -> list[BooleanSpace]:
    """
    Compute an optimized list of candidate states that is guaranteed
    to cover every attractor in the specified SD node.

    It should hold that for every attractor in the node's sub-space outside
    of the child sub-spaces (if known), the list contains at least one state
    from this attractor. However, it is not guaranteed that each candidate
    covers some attractor, so a non-empty list of candidates can
    still correspond to an empty list of attractors.

    The method starts by computing an NFVS of the network percolated to the
    node's space, and then assigns each NFVS node a value towards which
    it tends to. This transformation satisfies that the fixed-points in the
    resulting network cover the attractors in the original sub-space.

    Subsequently, the method can further optimize this set of candidates.

    Greedy ASP minification
    -----------------------
    Depending on which values are chosen for the NFVS nodes, a different set of
    candidates is obtained. This optimization method flips the initial NFVS node
    values and retains modifications that lead to a smaller candidate set, until
    a local minimum is achieved.

    Simulation minification
    -----------------------
    It is often possible to reduce the candidate set using a pseudo-random walk
    which eventually leads to one of the child spaces, or to one of the existing
    candidate states.

    Pint minification
    -----------------
    The tool `pint` can perform (incomplete) reachability check for any network
    state, hence we can also use it to reduce the list of candidate states.
    This option is only available if pint is actually installed.

    Parameters
    ----------
    sd : SuccessionDiagram
        The succession diagram in question.
    node_id : int
        An ID of an SD node for which we are computing the candidates.
    greedy_asp_minification: bool
        Whether to enable the iterative ASP-based minification.
    simulation_minification: bool
        Whether to enable simulation minification.
    pint_minification: bool
        Whether to enable pint minification.

    Returns
    -------
    list[BooleanSpace]
        The list of attractor candidate states.
    """

    if sd.config["debug"]:
        print(f"[{node_id}] Start computing attractor candidates.")

    node_data = sd.node_data(node_id)

    node_space = node_data["space"]

    if len(node_space) == sd.network.variable_count():
        if sd.config["debug"]:
            print(f"[{node_id}] > Attractor candidates done: node is a fixed-point.")
        return [node_space]

    node_nfvs = sd.node_percolated_nfvs(node_id, compute=True)

    if sd.config["debug"]:
        root_nfvs = sd.node_percolated_nfvs(sd.root(), compute=True)
        print(
            f"[{node_id}] > Percolated node NFVS contains {len(node_nfvs)} nodes (vs {len(root_nfvs)} in the root)."
        )

    # Compute the list of child stable-motifs if the node is expanded. Otherwise
    # "pretend" that there are no children.
    #
    # (We can safely ignore anything that is in these stable motifs, not just the
    # percolated child spaces, because it is either in the child space, or it percolates
    # to the child space, so it does not contain attractors)
    child_motifs_reduced = []
    if node_data["expanded"]:
        children = sd.node_successors(node_id, compute=False)
        child_motifs_reduced = [
            sd.edge_stable_motif(node_id, s, reduced=True) for s in children
        ]

    if node_data["skipped"]:
        # For skip nodes, it does not hold that the successors are the maximal subspaces.
        # This means that a skip node can intersect with some other SD node and that intersection
        # is not a subset of one of its children. In such case, we can use this intersection
        # to further simplify the attractor detection process.
        total_skip_nodes_applied = 0
        for n in sd.node_ids():
            n_data = sd.node_data(n)
            if is_subspace(node_space, n_data["space"]):
                # This means that (in a fully expanded SD), `node` would be a (transitive)
                # successor of `n`, which means that a result for `n` is "attractors in n
                # that are not in node". Hence, we can't use it to reason about `node`.
                # This is not a problem if the intersection of the two nodes is non-triviall,
                # because that means they have common successors (and those should be
                # solved separately), but the nodes themselves do not depend on each other.
                continue
            if n_data["attractor_candidates"] == [] or n_data["attractor_seeds"] == []:
                # This will create a lot of duplicates, but it seems to be better than
                # not doing it at all.
                common_subspace = intersect(node_space, n_data["space"])
                if common_subspace is not None:
                    reduced_subspace: BooleanSpace = {
                        k: v for k, v in common_subspace.items() if k not in node_space
                    }
                    child_motifs_reduced.append(reduced_subspace)
                    total_skip_nodes_applied += 1
        if sd.config["debug"]:
            print(
                f"[{node_id}] Extended child motifs with {total_skip_nodes_applied} skip-node intersections."
            )

    # Indicates that this space is either minimal, or has no computed successors.
    # In either case, the space must contain at least one attractor.
    node_is_pseudo_minimal = len(child_motifs_reduced) == 0

    if sd.config["debug"] and node_is_pseudo_minimal:
        print(
            f"[{node_id}] > The node has no children (i.e. it is minimal or unexpanded)."
        )

    if len(node_nfvs) == 0:
        # If the NFVS is empty, it means this node space has no complex attractors.
        # In a fully expanded SD, all fixed-points should be handled by the condition above,
        # so only non-minimal spaces remain and these can be resolved as empty.
        #
        # However, if SD is not fully expanded, it can contain further "pseudo-minimal"
        # spaces that would contain minimal traps if expanded. For those, there should
        # only be fixed-point attractors, but we need to identify
        # them, so we need to continue.

        # This tests the assumption that all fully expanded minimal
        # traps are handled before, and the only trap spaces with
        # an empty NFVS are unexpanded pseudo-minimal,
        # or expanded non-minimal.
        assert not sd.node_is_minimal(node_id)

        if not node_is_pseudo_minimal:
            if sd.config["debug"]:
                print(
                    f"[{node_id}] > Attractor candidates done: empty NFVS in a non-minimal space."
                )
            return []

    # First, we create a retained set and corresponding candidate set that is valid in the
    # percolated Boolean network. This process is *usually* quite fast, but
    # sometimes the candidate set explodes and is very large (even too large to enumerate).
    # In that case, we can try to simplify it greedily using the ASP solver (assuming this
    # is enabled).
    #
    # Overall, our strategy is the following:
    # 1. Compute a heuristic retained set based on the original NFVS method.
    # 2. If the ASP optimization is disabled, we generate the candidates from this retained
    #    set in full.
    # 3. If ASP optimization is enabled, we first test if there are at least 100 candidates
    #    for the initial retained set.
    #    3a. If there are fewer than 100 candidates, we just run one optimization pass on
    #        the initial retained set and report the result.
    #    3b. If there are more than 100 candidates, we probably need to optimize the retained
    #        set even further, in which case we start to recursively generate a new retained
    #        set that is optimized in each step to hopefully get a better result.

    pn_reduced = sd.node_percolated_petri_net(node_id, compute=True)
    bn_reduced = sd.node_percolated_network(node_id, compute=True)
    graph_reduced = AsynchronousGraph(bn_reduced)

    retained_set = make_heuristic_retained_set(
        graph_reduced, node_nfvs, child_motifs_reduced
    )

    if len(retained_set) == sd.network.variable_count() and node_is_pseudo_minimal:
        # If the retained set describes a fixed point, then only one attractor
        # is present in this space and it must contain the state described by the retained set.
        if sd.config["debug"]:
            print(
                f"[{node_id}] > Singular attractor found through fixed-point retained set. Done."
            )
        return [retained_set | node_space]

    if not greedy_asp_minification:
        candidate_states = compute_fixed_point_reduced_STG(
            pn_reduced,
            retained_set,
            avoid_subspaces=child_motifs_reduced,
            solution_limit=sd.config["attractor_candidates_limit"],
        )
        if len(candidate_states) == sd.config["attractor_candidates_limit"]:
            raise RuntimeError(
                f"Exceeded the maximum amount of attractor candidates ({sd.config['attractor_candidates_limit']}; see `SuccessionDiagramConfiguation.attractor_candidates_limit`)."
            )
        if sd.config["debug"]:
            print(
                f"[{node_id}] Computed {len(candidate_states)} candidate states without retained set optimization."
            )
    else:
        candidate_states = compute_fixed_point_reduced_STG(
            pn_reduced,
            retained_set,
            avoid_subspaces=child_motifs_reduced,
            solution_limit=sd.config["retained_set_optimization_threshold"],
        )

        if len(candidate_states) < sd.config["retained_set_optimization_threshold"]:
            # The candidate set is small. No need to overthink it. The original
            # retained set is probably quite good. However, we can still try to
            # optimize it further if it makes sense.
            if len(candidate_states) > 1 or (
                not node_is_pseudo_minimal and len(candidate_states) > 0
            ):
                if sd.config["debug"]:
                    print(
                        f"[{node_id}] Initial retained set generated {len(candidate_states)} candidates. Optimizing..."
                    )
                optimized = asp_greedy_retained_set_optimization(
                    sd,
                    node_id,
                    petri_net=pn_reduced,
                    retained_set=retained_set,
                    candidate_states=candidate_states,
                    avoid_dnf=child_motifs_reduced,
                )
                retained_set = optimized[0]
                candidate_states = optimized[1]
        else:
            # There seem to be many candidates, in which case it might be better
            # to optimize the retained set dynamically.
            if sd.config["debug"]:
                print(
                    f"[{node_id}] Initial retained set generated >{sd.config['retained_set_optimization_threshold']} candidates. Regenerate retained set."
                )
            retained_set = {}
            candidate_states = []
            for var in node_nfvs:
                retained_set[var] = 0
                candidate_states_zero = compute_fixed_point_reduced_STG(
                    pn_reduced,
                    retained_set,
                    avoid_subspaces=child_motifs_reduced,
                    solution_limit=sd.config["attractor_candidates_limit"],
                )

                if len(candidate_states_zero) <= len(candidate_states):
                    if sd.config["debug"]:
                        print(
                            f"[{node_id}] Chosen {var}=0 without increasing candidate count ({len(candidate_states_zero)}). {len(retained_set)}/{len(node_nfvs)} variables chosen."
                        )
                    candidate_states = candidate_states_zero
                    continue

                retained_set[var] = 1
                candidate_states_one = compute_fixed_point_reduced_STG(
                    pn_reduced,
                    retained_set,
                    avoid_subspaces=child_motifs_reduced,
                    solution_limit=len(candidate_states_zero),
                )

                if (
                    len(candidate_states_zero)
                    == sd.config["attractor_candidates_limit"]
                    and len(candidate_states_one)
                    == sd.config["attractor_candidates_limit"]
                ):
                    raise RuntimeError(
                        f"Exceeded the maximum amount of attractor candidates ({sd.config['attractor_candidates_limit']}; see `SuccessionDiagramConfiguation.attractor_candidates_limit`)."
                    )

                if len(candidate_states_one) <= len(candidate_states):
                    if sd.config["debug"]:
                        print(
                            f"[{node_id}] Chosen {var}=1 without increasing candidate count ({len(candidate_states_one)}). {len(retained_set)}/{len(node_nfvs)} variables chosen."
                        )
                    candidate_states = candidate_states_one
                    continue

                if len(candidate_states_zero) < len(candidate_states_one):
                    if sd.config["debug"]:
                        print(
                            f"[{node_id}] Chosen {var}=0 with better candidate count ({len(candidate_states_zero)}). {len(retained_set)}/{len(node_nfvs)} variables chosen."
                        )
                    candidate_states = candidate_states_zero
                    retained_set[var] = 0
                else:
                    if sd.config["debug"]:
                        print(
                            f"[{node_id}] Chosen {var}=1 with better candidate count ({len(candidate_states_one)}). {len(retained_set)}/{len(node_nfvs)} variables chosen."
                        )
                    candidate_states = candidate_states_one
                    retained_set[var] = 1

                # At this point, we know the candidate count increased and so we should
                # try to bring it back down.
                if (
                    len(candidate_states)
                    > sd.config["retained_set_optimization_threshold"]
                ):
                    if sd.config["debug"]:
                        print(f"[{node_id}] Optimizing partial retained set...")
                    optimized = asp_greedy_retained_set_optimization(
                        sd,
                        node_id,
                        petri_net=pn_reduced,
                        retained_set=retained_set,
                        candidate_states=candidate_states,
                        avoid_dnf=child_motifs_reduced,
                    )
                    retained_set = optimized[0]
                    candidate_states = optimized[1]

    # Terminate if done.
    if len(candidate_states) == 0:
        if sd.config["debug"]:
            print(f"[{node_id}] > Initial candidate set empty. Done.")
        return []
    if node_is_pseudo_minimal and len(candidate_states) == 1:
        if sd.config["debug"]:
            print(
                f"[{node_id}] > Single candidate found in (pseudo) minimal trap space. Done."
            )
        return [candidate_states[0] | node_space]

    if sd.config["debug"]:
        print(
            f"[{node_id}] > Attractor candidates from retained set: {len(candidate_states)}."
        )

    if simulation_minification:
        if sd.config["debug"]:
            print(f"[{node_id}] Start simulation minification...")
        avoid_children_symbolic = state_list_to_bdd(
            graph_reduced.symbolic_context(), child_motifs_reduced
        )

        # Here, we gradually increase the iteration count while
        # the candidate set is being actively reduced. If the simulation
        # cannot reduce any further states and exceeds the proposed budget,
        # we are done.
        iterations = 2**10
        max_budget = (
            sd.config["minimum_simulation_budget"] * bn_reduced.variable_count()
        )
        while len(candidate_states) > 0:
            if sd.config["debug"]:
                print(
                    f"[{node_id}] > Start simulation with {len(candidate_states)} states and simulation limit {iterations}."
                )
            reduced = run_simulation_minification(
                sd,
                node_id,
                graph_reduced,
                candidate_states,
                avoid_children_symbolic,
                max_iterations=iterations,
                simulation_seed=123,
            )

            if (
                len(reduced) == len(candidate_states)
                and (iterations * len(candidate_states)) > max_budget
            ):
                candidate_states = reduced
                break

            iterations = 2 * iterations
            candidate_states = reduced

            if len(candidate_states) == 1 and avoid_children_symbolic.is_false():
                break

        if sd.config["debug"]:
            print(f"[{node_id}] > Candidates after simulation: {len(candidate_states)}")

    # Terminate if done.
    if len(candidate_states) == 0:
        if sd.config["debug"]:
            print(f"[{node_id}] > Candidate set empty. Done.")
        return []
    if node_is_pseudo_minimal and len(candidate_states) == 1:
        if sd.config["debug"]:
            print(
                f"[{node_id}] > Single candidate found in (pseudo) minimal trap space. Done."
            )
        return [candidate_states[0] | node_space]

    if pint_minification and not pint_available:
        print("WARNING: Using `pint`, but `pint` is not installed. Skipping.")
    elif pint_minification and pint_available:
        if sd.config["debug"]:
            print(f"[{node_id}] Start `pint` minification...")

        children_bdd = state_list_to_bdd(
            graph_reduced.symbolic_context(), child_motifs_reduced
        )
        candidates_bdd = state_list_to_bdd(
            graph_reduced.symbolic_context(), candidate_states
        )
        avoid_bdd = children_bdd.l_or(candidates_bdd)

        filtered_states: list[BooleanSpace] = []
        for i, state in enumerate(candidate_states):
            state_bdd = state_to_bdd(graph_reduced.symbolic_context(), state)

            avoid_bdd = avoid_bdd.l_and_not(state_bdd)

            keep = True
            try:
                if biobalm._pint_reachability.pint_reachability(
                    pn_reduced, state, avoid_bdd, sd.config
                ):
                    keep = False
            except RuntimeError as e:
                assert str(e) == "Cannot verify."

            if keep:
                avoid_bdd = avoid_bdd.l_or(state_bdd)
                filtered_states.append(state)

            if sd.config["debug"]:
                print(
                    f"[{node_id}] > `pint` {i + 1}/{len(candidate_states)}: eliminated: {not keep}, retained: {len(filtered_states)}."
                )

        candidate_states = filtered_states

        if sd.config["debug"]:
            print(f"[{node_id}] > Candidates after `pint`: {len(candidate_states)}")

    # Finally, we need to augment each candidate state with the
    # fixed values from the node space, because until now, we only considered
    # the reduced/percolated network.
    candidate_states = [x | node_space for x in candidate_states]

    return candidate_states


def run_simulation_minification(
    sd: SuccessionDiagram,
    node_id: int,
    graph: AsynchronousGraph,
    candidate_states: list[BooleanSpace],
    avoid_bdd: Bdd,
    max_iterations: int,
    simulation_seed: int,
) -> list[BooleanSpace]:
    """
    A fast but incomplete method for eliminating spurious attractor candidates
    based on stochastic simulation.

    Parameters
    ----------
    node_id : int
        The ID of the associated SD node. This is only for logging progress.
    graph : biodivine_aeon.AsynchronousGraph
        The symbolic encoding of the *percolated* network dynamics.
    candidate_states: list[BooleanSpace]
        The list of currently considered candidates that is to be reduced.
    avoid_bdd: biodivine_aeon.Bdd
        A symbolic representation of the states/spaces that are to be ignored.
        If any of these states is reachable by a candidate, that candidate is
        safe to ignore as well.
    max_iterations: int
        The number of steps performed by the simulator.
    simulation_seed: int
        The seed value for the random walk simulator.

    Returns
    -------
    list[BooleanSpace]
        The optimized list of candidate states.
    """

    # A random generator initialized with a fixed seed. Ensures simulation
    # is randomized but deterministic.
    generator = random.Random(simulation_seed)

    # Retrieve the symbolic variables that encode each network variable and
    # associate them with the corresponding network's update functions.
    symbolic_ctx = graph.symbolic_context()
    network_vars = graph.network_variables()
    symbolic_vars: list[BddVariable] = []
    for var in network_vars:
        s_var = symbolic_ctx.find_network_bdd_variable(var)
        assert s_var is not None
        symbolic_vars.append(s_var)
    update_functions = {
        symbolic_vars[i]: graph.mk_update_function(network_vars[i])
        for i in range(len(network_vars))
    }

    # Use stochastic simulation to prune the set of candidate states.
    # We use different simulation approach depending on whether this space
    # is a minimal trap or not.
    if not avoid_bdd.is_false():
        candidates_bdd = state_list_to_bdd(symbolic_ctx, candidate_states)
        filtered_candidates: list[BooleanSpace] = []

        for i, state in enumerate(candidate_states):
            if i % 100 == 99 and sd.config["debug"]:
                print(
                    f"[{node_id}] > Simulation progress: {i + 1}/{len(candidate_states)}"
                )

            # Remove the state from the candidates. If we can prove that is
            # is not an attractor, we will put it back.
            state_bdd = graph.mk_subspace(state).to_bdd()
            candidates_bdd = candidates_bdd.l_and_not(state_bdd)

            # Convert the state BDD to a valuation object. This is like
            # a dictionary, but updates much faster because it is managed
            # by rust and indexed by symbolic variable IDs and not names.
            #
            # (The `state_bdd` should be a singleton, so it's ok to just
            # take the first valuation)
            simulation = state_bdd.valuation_first()

            is_valid_candidate = True
            for _ in range(max_iterations):
                # Advance all variables by one step in random order.
                generator.shuffle(symbolic_vars)
                for s_var in symbolic_vars:
                    # We can just "call" the BDD as a normal function with
                    # a valuation as the input.
                    step = update_functions[s_var](simulation)
                    simulation[s_var] = step

                # To check for set inclusion, we can also just "eval" the
                # BDD, since it evaluates to true for every state that is
                # in the encoded set.

                if candidates_bdd(simulation):
                    # The state can reach some other state in the candidate
                    # set. This does not mean it cannot be an attractor, but
                    # it means it is sufficient to keep considering
                    # the remaining candidates.
                    is_valid_candidate = False
                    break

                if avoid_bdd(simulation):
                    # The state can reach some other state in the avoid set,
                    # which means it cannot be an attractor in this subspace.
                    is_valid_candidate = False
                    break

            if is_valid_candidate:
                # If we cannot rule out the candidate, we have to put it back
                # into the candidate set.

                # However, we want to use the new state (after simulation),
                # because this one is probably "closer" to the goals that we
                # want to eventually get to.

                candidates_bdd = candidates_bdd.l_or(Bdd(simulation))
                filtered_candidates.append(valuation_to_state(symbolic_ctx, simulation))

        return filtered_candidates
    else:
        # If the avoid set is empty, it means that this is a pseudo-minimal space
        # and hence we (a) don't have to check if we reached the avoid set,
        # and (b) we can stop with one candidate instead of zero.

        candidates_bdd = state_list_to_bdd(symbolic_ctx, candidate_states)
        printed: set[int] = set()
        for i in range(max_iterations):
            progress = int((i * len(candidate_states)) / max_iterations)
            if (
                progress % 100 == 99
                and sd.config["debug"]
                and (progress not in printed)
            ):
                printed.add(progress)
                print(
                    f"[{node_id}] > Simulation progress: {progress + 1}/{len(candidate_states)}"
                )

            generator.shuffle(symbolic_vars)
            new_candidates_bdd = symbolic_ctx.mk_constant(False)

            # Goes through all the states in the candidates BDD, updates
            # them, and puts them in the new BDD assuming they are not there yet.
            for state_val in candidates_bdd.valuation_iterator():
                candidates_bdd = candidates_bdd.l_and_not(Bdd(state_val))

                simulation = state_val
                for s_var in symbolic_vars:
                    step = update_functions[s_var](simulation)
                    simulation[s_var] = step

                if candidates_bdd(simulation) or new_candidates_bdd(simulation):
                    # We have reached one of the other candidates.
                    continue

                new_candidates_bdd = new_candidates_bdd.l_or(Bdd(simulation))

            candidates_bdd = new_candidates_bdd

            if candidates_bdd.cardinality() <= 1:
                break

        # Finally, we convert the candidate states back into something that
        # biobalm understands.
        filtered_candidates = []
        for state_val in candidates_bdd.valuation_iterator():
            filtered_candidates.append(valuation_to_state(symbolic_ctx, state_val))

        return filtered_candidates


def asp_greedy_retained_set_optimization(
    sd: SuccessionDiagram,
    node_id: int,
    petri_net: DiGraph,
    retained_set: BooleanSpace,
    candidate_states: list[BooleanSpace],
    avoid_dnf: list[BooleanSpace],
) -> tuple[BooleanSpace, list[BooleanSpace]]:
    """
    Takes a Boolean network encoded as a Petri net and a candidate retained set.
    Then, performs a greedy optimization process leading to a locally optimal
    version of the original retained set that results in the least
    amount of candidate states.


    Parameters
    ----------
    node_id : int
        The ID of the associated SD node. This is only for logging progress.
    petri_net : networkx.DiGraph
        The Petri net encoding of the *percolated* network dynamics.
    retained_set : BooleanSpace
        The retained set that will be the initial point of the optimization.
    candidate_states: list[BooleanSpace]
        The list of candidates that are valid for the given `retained_set`, so
        that they don't need to be recomputed.
    avoid_dnf: list[BooleanSpace]
        The list of subspaces in the given network in which candidate
        states can be ignored.

    Returns
    -------
    tuple[BooleanSpace, list[BooleanSpace]]
        The optimized reatined set, together with the list of candidates that
        are valid for this retained set.
    """
    done = False
    while not done:
        done = True
        for var in retained_set:
            # Standrad termination checks.
            if len(candidate_states) == 0:
                return (retained_set, [])
            if len(avoid_dnf) == 0 and len(candidate_states) == 1:
                # If we are not avoiding anything, this is a (pseudo)
                # minimal space and a single candidate is good enough.
                return (retained_set, candidate_states)

            # Flip the value of the specified variable to see if the resulting
            # candidate set is smaller.
            retained_set_2 = retained_set.copy()
            retained_set_2[var] = cast(Literal[0, 1], 1 - retained_set_2[var])
            candidate_states_2 = compute_fixed_point_reduced_STG(
                petri_net,
                retained_set_2,
                avoid_subspaces=avoid_dnf,
                # We don't need all solutions if the result isn't smaller.
                solution_limit=len(candidate_states),
            )
            if len(candidate_states_2) < len(candidate_states):
                retained_set = retained_set_2
                candidate_states = candidate_states_2
                done = False
                if sd.config["debug"]:
                    print(
                        f"[{node_id}] > Candidate states optimized to {len(candidate_states)}."
                    )
    return (retained_set, candidate_states)


def make_heuristic_retained_set(
    graph: AsynchronousGraph, nfvs: list[str], avoid_dnf: list[BooleanSpace]
) -> BooleanSpace:
    """
    Calculate the retained set for a given `space` based on heuristic criteria.

    The retained set is technically a space-like object that describes the
    variables which have to be fixed in order for the network to lose all
    complex attractors. However, note that this really means changing the update
    functions of these variables. This is not a trap space that only contains
    fixed-points, but a description of how the network must be modified to
    eliminate all complex attractors.

    The construction guarantees that any complex attractor of the old
    network will manifest as at least one fixed-point in the new network. But this
    network modification is not implemented here. This method only generates
    the necessary list of variables and values.

    Finally, note that the method uses a heuristic to select values
    that should ideally lead to the least amount of fixed-points in
    the modified network, but there can be other, more optimal retained sets.

    Parameters
    ----------
    graph : AsynchronousGraph
        The symbolic update functions stored as an `AsynchronousGraph` object
        from the `biodivine_aeon` library.
    nfvs : list[str]
        The list of variables that represent an NFVS of the network encoded
        by `graph`.
    avoid_dnf : list[BooleanSpace]
        A list of :class:`BooleanSpace<biobalm.types.BooleanSpace>` objects
        describing the child spaces of the specified network. Only attractors that
        are not in these child spaces are considered. If no child spaces are
        provided, then all attractors are considered.

    Returns
    -------
    BooleanSpace
        A :class:`BooleanSpace<biobalm.types.BooleanSpace>` object describing the
        retained set.
    """
    retained_set: BooleanSpace = {}

    # First, if there are any child spaces present, we extend the retained set
    # with the values from the one that has the least amount of fixed variables
    # shared with the NFVS.
    if len(avoid_dnf) > 0:
        # Find the child space that has the fewest nodes in common with the NFVS:
        least_common_child_space = avoid_dnf[0]
        least_common_nodes = len(set(least_common_child_space) & set(nfvs))
        for child_space in avoid_dnf:
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
