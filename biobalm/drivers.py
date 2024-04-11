"""
A module which provides functions for processing the stable motif driver nodes of a Boolean network.
"""

from __future__ import annotations

from biodivine_aeon import AsynchronousGraph, BooleanNetwork

from biobalm.space_utils import percolate_space_strict
from biobalm.types import BooleanSpace


def find_single_node_LDOIs(
    network: AsynchronousGraph | BooleanNetwork,
) -> dict[tuple[str, int], BooleanSpace]:
    """
    Finds the the logical domain of influence for every single node state.

    The logical domain of influence is equivalent to the set of node states that
    become fixed after percolation. See
    :func:`percolate_space<biobalm.space_utils.percolate_space>` for details.

    This operation requires the symbolic update functions provided by an
    `AsynchronousGraph`. If you provide a `BooleanNetwork`, the
    `AsynchronousGraph` will be created automatically, but this can incur
    additional overhead.

    Parameters
    ----------
    network : AsynchronousGraph | BooleanNetwork
        The symbolic update functions stored as a `AsynchronousGraph` or
        `BooleanNetwork` object from the `biodivine_aeon` library.

    Returns
    -------
    dict[tuple[str, int], BooleanSpace]
        The LDOI for every one-node state. The keys are `(variable, value)`
        tuples, where `variable` is the name of the node that is fixed in the
        state `value` (`0` or `1`). The keys are
        :class:`BooleanSpace<biobalm.types.BooleanSpace>` objects, which are dictionaries
        of node values describing the node states that become fixed as a result
        of percolation.

    Example
    _______
    >>> import biobalm
    >>> sd = biobalm.SuccessionDiagram.from_rules('A, A\\nB, A')
    >>> ldois = biobalm.drivers.find_single_node_LDOIs(sd.network)
    >>> for k in sorted(ldois):
    ...     print(f'{k} ==> {sorted(ldois[k].items())}')
    ...
    ('A', 0) ==> [('A', 0), ('B', 0)]
    ('A', 1) ==> [('A', 1), ('B', 1)]
    ('B', 0) ==> []
    ('B', 1) ==> []
    """
    if isinstance(network, BooleanNetwork):
        network = AsynchronousGraph(network)

    LDOIs: dict[tuple[str, int], BooleanSpace] = {}
    for var in network.network_variable_names():
        fn_bdd = network.mk_update_function(var)
        if fn_bdd.is_true() or fn_bdd.is_false():
            # Skip constant nodes.
            continue

        LDOIs[(var, 0)] = percolate_space_strict(network, {var: 0})
        LDOIs[(var, 1)] = percolate_space_strict(network, {var: 1})

    return LDOIs


def find_single_drivers(
    target_subspace: BooleanSpace,
    network: AsynchronousGraph | BooleanNetwork,
    LDOIs: dict[tuple[str, int], BooleanSpace] | None = None,
) -> set[tuple[str, int]]:
    """
    Find all one-node drivers for a target subspace.

    This operation requires the symbolic update functions provided by an
    `AsynchronousGraph`. If you provide a `BooleanNetwork`, the
    `AsynchronousGraph` will be created automatically, but this can incur
    additional overhead.

    Parameters
    ----------
    target_subspace : BooleanSpace
        A :class:`BooleanSpace<biobalm.types.BooleanSpace>` object describing the
        target trap space.
    network : AsynchronousGraph | BooleanNetwork
       The symbolic update functions stored as a `AsynchronousGraph` or
        `BooleanNetwork` object from the `biodivine_aeon` library.
    LDOIs : dict[tuple[str, int], BooleanSpace] | None, optional
        A dictionary encoding the LDOI for every one-node state to be
        considered. The keys should be `(variable, value)` tuples, where
        `variable` is the name of the node that is fixed in the state `value`
        (`0` or `1`). The keys are
        :class:`BooleanSpace<biobalm.types.BooleanSpace>` objects, which are
        dictionaries of node values describing the node states that become fixed
        as a result of percolation. If not provided, the LDOI will be computed
        automatically for each node state using :func:`find_single_node_LDOIs`.


    Returns
    -------
    set[tuple[str, int]]
        A set of driver node states that result in the target subspace, with
        each node state represented as a node name and value (`0` or `1`).

    Example
    -------
    >>> import biobalm
    >>> sd = biobalm.SuccessionDiagram.from_rules('A, A\\nB, A')
    >>> drivers = biobalm.drivers.find_single_drivers({'B': 0}, sd.network)
    >>> sorted(list(drivers))
    [('A', 0), ('B', 0)]
    """
    if isinstance(network, BooleanNetwork):
        network = AsynchronousGraph(network)

    if LDOIs is None:
        LDOIs = find_single_node_LDOIs(network)

    drivers: set[tuple[str, int]] = set()
    for fix, LDOI in LDOIs.items():
        if target_subspace.items() <= (LDOI.items() | {fix}):
            drivers.add(fix)

    return drivers
