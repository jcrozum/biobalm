from biodivine_aeon import (
    BooleanNetwork,
    FixedPoints,
    AsynchronousGraph,
)

from balm.petri_net_translation import network_to_petrinet
from balm.trappist_core import compute_fixed_point_reduced_STG, trappist
from balm.types import BooleanSpace


def remove_static_constraints(network: BooleanNetwork) -> BooleanNetwork:
    """
    A method that removes all information about regulation monotonicity and
    essentiality from the given `BooleanNetwork`.

    This is mostly done to allow handling of randomly generated or otherwise
    machine pre-processed files that can contain subtle logical redundancies
    that AEON would otherwise detect as warnings.
    """
    bn = BooleanNetwork(network.variable_names())
    for reg in network.regulations():
        reg['essential'] = False
        reg['sign'] = None
        bn.add_regulation(reg)
    
    for var in network.variables():
        bn.set_update_function(var, network.get_update_function(var))

    return bn


def test_network_minimum_traps(network_file: str):
    bn = remove_static_constraints(BooleanNetwork.from_file(network_file))
    stg = AsynchronousGraph(bn)

    min_max_traps = trappist(bn, problem="min") + trappist(bn, problem="max")

    pn = network_to_petrinet(bn)
    min_max_traps_pre_encoded = trappist(pn, problem="min") + trappist(
        pn, problem="max"
    )

    assert min_max_traps == min_max_traps_pre_encoded

    # We have no way of knowing if a trap is minimal/maximal, but we can still
    # verify that it is indeed a trap.
    for trap in min_max_traps:
        # Then a proper symbolic check that should be reliable every time.
        symbolic_space = stg.mk_subspace(trap)
        if stg.post(symbolic_space).is_subset(symbolic_space):
            continue
        raise Exception(f"Failed on {network_file}: {trap} is not a trap space.")


def test_network_fixed_points(network_file: str):
    # Verify that the fixed-points of the test models are the same
    # as when computing using BDDs.
    bn = remove_static_constraints(BooleanNetwork.from_file(network_file))
    stg = AsynchronousGraph(bn)

    symbolic_fixed_points = FixedPoints.symbolic(stg, stg.mk_unit_colored_vertices())
    trappist_fixed_points = trappist(bn, problem="fix")
    for fixed_point in trappist_fixed_points:
        # Convert trappist result to a symbolic singleton set.
        vertex = stg.mk_subspace(fixed_point)
        # Check that the fixed-point has been found, and remove it.
        assert vertex.is_subset(
            symbolic_fixed_points
        ), f"Failed on {network_file}: {fixed_point} is not in symbolic fixed points."
        symbolic_fixed_points = symbolic_fixed_points.minus(vertex)
    # In the end, all fixed-points must have been found.
    assert (
        symbolic_fixed_points.is_empty()
    ), f"Failed on {network_file}: Some symbolic fixed points not detected by trappist."


def test_network_fixed_point_reduced_STG():
    # Validate the function for computing fixed points of the reduced STG
    # on a single small input.

    bn = BooleanNetwork.from_bnet(
        """
        x1, (x1 & x2) | (!x1 & !x2)
        x2, (x1 & x2) | (!x1 & !x2)
    """
    )

    petri_net = network_to_petrinet(bn)

    avoid_subspace_1: BooleanSpace = {"x1": 1, "x2": 1}
    avoid_subspace_2: BooleanSpace = {}
    avoid_subspace_3: BooleanSpace = {"x2": 1}

    ensure_subspace_1: BooleanSpace = {}
    ensure_subspace_2: BooleanSpace = {"x1": 0, "x2": 0}

    retained_set: BooleanSpace = {"x1": 0, "x2": 0}
    candidate_set = compute_fixed_point_reduced_STG(petri_net, retained_set)
    assert len(candidate_set) == 2  # candidate_set = {00, 11}

    candidate_set = compute_fixed_point_reduced_STG(
        petri_net,
        retained_set,
        ensure_subspace=ensure_subspace_1,
        avoid_subspaces=[avoid_subspace_1],
    )
    assert len(candidate_set) == 1  # candidate_set = {00}

    candidate_set = compute_fixed_point_reduced_STG(
        petri_net,
        retained_set,
        ensure_subspace=ensure_subspace_1,
        avoid_subspaces=[avoid_subspace_2],
    )
    assert len(candidate_set) == 0  # candidate_set = empty

    candidate_set = compute_fixed_point_reduced_STG(
        petri_net, retained_set, ensure_subspace=ensure_subspace_2
    )
    assert len(candidate_set) == 1  # candidate_set = {00}

    retained_set = {"x1": 1, "x2": 1}
    candidate_set = compute_fixed_point_reduced_STG(petri_net, retained_set)
    assert len(candidate_set) == 3  # candidate_set = {01, 10, 11}

    candidate_set = compute_fixed_point_reduced_STG(
        petri_net,
        retained_set,
        ensure_subspace=ensure_subspace_1,
        avoid_subspaces=[avoid_subspace_1],
    )
    assert len(candidate_set) == 2  # candidate_set = {01, 10}

    candidate_set = compute_fixed_point_reduced_STG(
        petri_net,
        retained_set,
        ensure_subspace=ensure_subspace_1,
        avoid_subspaces=[avoid_subspace_2],
    )
    assert len(candidate_set) == 0  # candidate_set = empty

    candidate_set = compute_fixed_point_reduced_STG(
        petri_net,
        retained_set,
        ensure_subspace=ensure_subspace_1,
        avoid_subspaces=[avoid_subspace_3],
    )
    assert len(candidate_set) == 1  # candidate_set = {10}
