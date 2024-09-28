"""
This module implementes methods for permanent target control of Boolean networks
based on the structure of a succession diagram.
"""

from __future__ import annotations

from functools import reduce
from itertools import combinations, product
from typing import Iterator, Literal, cast

import networkx as nx  # type: ignore
from biodivine_aeon import AsynchronousGraph, BooleanNetwork

from biobalm.space_utils import intersect, is_subspace, percolate_space
from biobalm.succession_diagram import SuccessionDiagram
from biobalm.types import BooleanSpace, ControlOverrides, SubspaceSuccession


class Intervention:
    def __init__(
        self,
        control: list[ControlOverrides],
        strategy: Literal["internal", "all"],
        succession: SubspaceSuccession,
    ):
        """A class for encoding an intervention to a network to reach a target subspace.

        Generally, this class is created by the
        :func:`succession_control<biobalm.control.succession_control>` function,
        which returns a list of `Intervention` objects. Manipulating the
        contents of these objects is only recommended for advanced use cases.
        Typically, it is sufficient to print this object to see a human-readable
        explanation of how to interpret the intervention.

        Note that two interventions are considered to be equal if they act on the same succession with
        equal controls applied to each subspace in the succession. Thus, two
        interventions that override the same nodes in the same way can be
        unequal, even when applied to the same subspaces. This is because each
        subspace imposes a duration condition on the overrides, i.e., the
        override must be maintained until the succession subspace is reached.
        Changing the order of the subspaces considered can alter this stop
        condition.

        Parameters
        ----------
        control : list[ControlOverrides]
            The :class:`ControlOverrides<biobalm.types.ControlOverrides>` objects,
            in order, that are applied. The order of the list corresponds the
            order of the subspaces in the succession. Each
            :class:`ControlOverrides<biobalm.types.ControlOverrides>` object
            represens a list of overrides, stored as a dictionary of node-value
            pairs, that drive the system to the corresponding subspace. Each
            :class:`ControlOverrides<biobalm.types.ControlOverrides>` object is
            sorted by key value (i.e., alphabetically) upon creation of the
            `Intervention` object to maintain a canonical ordering.
        strategy : str
            Either "internal" or "all"; "internal" means that the
            :class:`ControlOverrides<biobalm.types.ControlOverrides>`
        succession : SubspaceSuccession
            A sequence of subspaces that are targeted by the corresponding
            entries of `control`.

        Example
        -------
        >>> import biobalm
        >>> sd = biobalm.SuccessionDiagram.from_rules(
        ...     \"\"\"
        ...     A, B & C
        ...     B, A & C
        ...     C, A & B
        ...     \"\"\"
        ...     )
        >>> target = {"A": 1, "B": 1, "C": 1}
        >>> interventions = biobalm.control.succession_control(sd, target)
        >>> intervention = interventions[0] # only one in this case
        >>> intervention.control
        [[{'A': 1, 'B': 1}, {'A': 1, 'C': 1}, {'B': 1, 'C': 1}]]
        >>> intervention.strategy
        'internal'
        >>> intervention.succession
        [{'A': 1, 'B': 1, 'C': 1}]
        >>> intervention.successful
        True
        >>> print(intervention)
        Intervention is SUCCESSFUL operating on
        {'A': 1, 'B': 1, 'C': 1}
        override
        ({'A': 1, 'B': 1} or {'A': 1, 'C': 1} or {'B': 1, 'C': 1})
        """

        # we store the control in a canonical (sorted) representation
        self._control: list[ControlOverrides] = []
        for c in control:
            cs = sorted(map(lambda x: sorted(x.items()), c))
            self._control.append(list(map(dict, cs)))  # type: ignore

        self._strategy = strategy
        self._succession = succession
        self._successful = not any(not c for c in control)

    @property
    def control(self):
        return self._control

    @property
    def strategy(self):
        return self._strategy

    @property
    def succession(self):
        return self._succession

    @property
    def successful(self):
        return self._successful

    def __eq__(self, other: object):
        if not isinstance(other, Intervention):
            return False

        if self.succession != other.succession:
            return False

        if len(self.control) != len(other.control):
            return False

        for d1, d2 in zip(self.control, other.control):
            if not controls_are_equal(d1, d2):
                return False

        return True

    def __repr__(self):
        return (
            f"Intervention("
            f"{self.control},"
            f"{self.strategy},"
            f"{self.succession},"
            f"{self.successful})"
        )

    def __str__(self):
        succession_string = (
            f"Intervention is {'' if self.successful else 'UN'}SUCCESSFUL operating on\n"
            + "\n".join(map(str, self.succession))
            + "\noverride\n"
        )
        if self.strategy == "internal":
            return succession_string + " and then\n".join(
                f"({' or '.join(map(str, motif_control))})"
                for motif_control in self.control
            )
        elif self.strategy == "all":
            return succession_string + " temporarily, and then\n".join(
                f"({' or '.join(map(str, motif_control))})"
                for motif_control in self.control
            )
        else:
            return "unknown strategy: " + self.__repr__()

    def all_control_strategies(self) -> Iterator[ControlOverrides]:
        """
        Returns all possible combinations of `ControlOverrides` sequences that
        can be used to execute this `Intervention`.

        Internally, an intervention consists of multiple control steps that
        need to be taken sequentially. For each step in the sequence, an intervention
        can have multiple options of how to execute it. With this method,
        we can generate the actual sequences that arise by combining all the
        available options for each step.
        """
        return map(lambda x: list(x), product(*self._control))


def succession_control(
    succession_diagram: SuccessionDiagram,
    target: BooleanSpace,
    strategy: Literal["internal", "all"] = "internal",
    max_drivers_per_succession_node: int | None = None,
    forbidden_drivers: set[str] | None = None,
    successful_only: bool = True,
    skip_feedforward_successions: bool = False,
) -> list[Intervention]:
    """
    Performs succession-diagram control to reach a target subspace.

    Parameters
    ----------
    succession_diagram : SuccessionDiagram
        The succession diagram from which successions and rules will be
        extracted.
    target : BooleanSpace
        The target subspace.
    strategy : str, optional
        The searching strategy to use to look for driver nodes. Options are
        'internal' (default) and 'all'.
    max_drivers_per_succession_node: int | None = None,
        The maximum number of drivers that will be tested for a succession
        diagram node. If `None`, then a number of drivers up to the size of the
        succession diagram node's stable motif will be tested
    forbidden_drivers: set[str] | None
        A set of forbidden drivers that will not be overridden for control. If
        `None`, then all nodes are candidates for control.
    successful_only: bool
        Whether to only return successful interventions (default: `True`).

    Returns
    -------
    list[Intervention]
        A list of control :class:`Intervention<biobalm.control.Intervention>`
        objects. Note that when `successful_only` is `False`, returned
        interventions may be unsuccessful if `max_drivers_per_succession_node`
        is set too small, or crucial nodes are included in `forbidden_drivers`.
        To test, examine the `successful` property of the intervention.

    Example
    -------
    >>> import biobalm
    >>> from biobalm.control import succession_control
    >>> sd = biobalm.SuccessionDiagram.from_rules(
    ...     \"\"\"
    ...     S, S
    ...     A, S | B
    ...     B, A
    ...     C, A | D
    ...     D, C
    ...     E, false
    ...     \"\"\"
    ...     )
    >>> target = {"S": 0, "E": 0, "A": 0, "B": 0, "C": 1, "D": 1}
    >>> interventions = succession_control(sd, target, forbidden_drivers = {"A"})
    >>> interventions.sort(key=lambda x: len(x.control)) # to maintain fixed order
    >>> for i in interventions:
    ...     print(f'{i}\\n'+'-'*20)
    ...
    Intervention is SUCCESSFUL operating on
    {'S': 0}
    {'A': 0, 'B': 0}
    {'C': 1, 'D': 1}
    override
    ({'S': 0}) and then
    ({'B': 0}) and then
    ({'C': 1} or {'D': 1})
    --------------------
    Intervention is SUCCESSFUL operating on
    {'S': 0}
    {'C': 1, 'D': 1}
    {'A': 0, 'B': 0}
    override
    ({'S': 0}) and then
    ({'C': 1} or {'D': 1}) and then
    ({'B': 0})
    --------------------

    Example
    -------
    >>> import biobalm
    >>> from biobalm.control import succession_control
    >>> sd = biobalm.SuccessionDiagram.from_rules(
    ...         \"\"\"
    ...     S, S
    ...     A, S | B
    ...     B, A
    ...     C, A | D
    ...     D, C
    ...     E, false
    ...     \"\"\"
    ...     )
    >>> target = {"S": 0, "E": 0, "A": 0, "B": 0, "C": 1, "D": 1}
    >>> interventions = succession_control(sd, target, strategy = "all")
    >>> interventions.sort(key=lambda x: len(x.control)) # to maintain fixed order
    >>> for i in interventions:
    ...     print(f'{i}\\n'+'-'*20)
    ...
    Intervention is SUCCESSFUL operating on
    {'S': 0}
    {'A': 0, 'B': 0}
    {'C': 1, 'D': 1}
    override
    ({'S': 0}) temporarily, and then
    ({'A': 0} or {'B': 0}) temporarily, and then
    ({'C': 1} or {'D': 1})
    --------------------
    Intervention is SUCCESSFUL operating on
    {'S': 0}
    {'C': 1, 'D': 1}
    {'A': 0, 'B': 0}
    override
    ({'S': 0}) temporarily, and then
    ({'A': 1} or {'B': 1} or {'C': 1} or {'D': 1}) temporarily, and then
    ({'A': 0} or {'B': 0})
    --------------------
    """
    interventions: list[Intervention] = []

    successions = successions_to_target(
        succession_diagram,
        target=target,
        expand_diagram=True,
        skip_feedforward_successions=skip_feedforward_successions,
    )

    for succession in successions:
        controls = drivers_of_succession(
            succession_diagram.symbolic,
            succession,
            strategy=strategy,
            max_drivers_per_succession_node=max_drivers_per_succession_node,
            forbidden_drivers=forbidden_drivers,
        )
        intervention = Intervention(controls, strategy, succession)

        if not successful_only or intervention.successful:
            interventions.append(intervention)

    return interventions


def successions_to_target(
    succession_diagram: SuccessionDiagram,
    target: BooleanSpace,
    expand_diagram: bool = True,
    skip_feedforward_successions: bool = False,
) -> list[SubspaceSuccession]:
    """Find lists of nested trap spaces (successions) that lead to the
    specified target subspace.

    Generally, it is not necessary to call this function directly, as it is
    automatically invoked by
    :func:`succession_control<biobalm.control.succession_control>`. It is
    primarily provided in the public API for testing and benchmarking purposes,
    or in the case that the user wants to implement a custom strategy to
    identify succession drivers rather than relying on
    :func:`drivers_of_succession<biobalm.control.drivers_of_succession>`.

    Parameters
    ----------
    succession_diagram : SuccessionDiagram
        The succession diagram from which successions will be extracted.
    target : BooleanSpace
        The target subspace.
    expand_diagram: bool
        Whether to ensure that the succession diagram is expanded enough to
        capture all paths to the target (default: True).
    skip_feedforward_successions: bool
        Whether to skip redundant successions (default: False). Skipping these
        can reduce the number of interventions to test, yielding performance
        improvements, but can also cause the algorithm to miss some valid
        interventions, particularly in cases when the order of intervention
        application is important.

    Returns
    -------
    list[SubspaceSuccession]
        A list of successions, where each succession is a list of sequentially
        nested trap spaces that specify the target.
    """
    successions: list[SubspaceSuccession] = []
    # Tracks the combined perturbation that needs to be applied for the whole
    # succession to take effect. We use this to detect which successions are
    # redundant when they can be replaced by a succession with a subset
    # signature. Used when skip_feedforward_successions is True
    succession_signatures: list[BooleanSpace] = []

    # expand the succession_diagram toward the target
    if expand_diagram:
        succession_diagram.expand_to_target(
            target=target,
        )

    # these contradict the target or have a motif avoidant attractor that isn't
    # in full agreement with the target
    hot_lava_nodes: set[int] = set()

    descendant_map: dict[int, set[int]] = {}
    for s in succession_diagram.node_ids():
        is_consistent = intersect(succession_diagram.node_data(s)["space"], target)
        is_goal = is_subspace(succession_diagram.node_data(s)["space"], target)
        is_minimal = succession_diagram.node_is_minimal(s)
        if not is_consistent or (not is_goal and is_minimal):
            hot_lava_nodes.add(s)

        descendant_map[s] = set(nx.descendants(succession_diagram.dag, s))  # type: ignore
        descendant_map[s].add(s)  # for our purposes, s is its own descendant
    found_valid_target_node = False
    for s in succession_diagram.node_ids():
        # a node is a valid end point if all
        # 1) its descendents (including itself) are not "hot lava"
        # 2) it has a parent that is or can reach "hot lava" (otherwise, we can
        # just control to the parent) note that all nodes are either cold lava
        # or hot lava, but not both or neither
        if descendant_map[s] & hot_lava_nodes:
            continue
        found_valid_target_node = True
        if not any(
            descendant_map[p] & hot_lava_nodes
            for p in succession_diagram.dag.predecessors(s)  # type: ignore
        ):
            continue

        for path in cast(
            list[list[int]],
            nx.all_simple_paths(  # type: ignore
                succession_diagram.dag,
                source=succession_diagram.root(),
                target=s,
            ),
        ):
            motif_list = [
                succession_diagram.edge_all_stable_motifs(x, y, reduced=True)
                for x, y in zip(path[:-1], path[1:])
            ]
            for succession_tuple in product(*motif_list):
                succession = list(succession_tuple)
                if skip_feedforward_successions:
                    signature = reduce(lambda x, y: x | y, succession)
                    # First, check if any existing successions can be eliminated
                    # because they are redundant w.r.t. to this succession.
                    # (`reversed` is important here, because that way a delete
                    # only impacts indices that we already processed)
                    skip_completely = False
                    for i in reversed(range(len(succession_signatures))):
                        existing_signature = succession_signatures[i]
                        if is_subspace(signature, existing_signature):
                            # The current `path` is already superseded by a path in successions.
                            skip_completely = True
                            break
                        if is_subspace(existing_signature, signature):
                            # A path in successions is made redundant by the current path.
                            del succession_signatures[i]
                            del successions[i]
                    if skip_completely:
                        continue

                    succession_signatures.append(signature)
                successions.append(succession)

    if found_valid_target_node and len(successions) == 0:
        successions = [[]]

    return successions


def drivers_of_succession(
    bn: BooleanNetwork | AsynchronousGraph,
    succession: list[BooleanSpace],
    strategy: str = "internal",
    max_drivers_per_succession_node: int | None = None,
    forbidden_drivers: set[str] | None = None,
) -> list[ControlOverrides]:
    """
    Find driver nodes of a list of sequentially nested trap spaces.

    Generally, it is not necessary to call this function directly, as it is
    automatically invoked by
    :func:`succession_control<biobalm.control.succession_control>`. It is primarily
    provided in the public API for testing and benchmarking purposes.

    Parameters
    ----------
    bn : BooleanNetwork | AsynchronousGraph
        The network to analyze, which contains the Boolean update functions. Ideally,
        the network should be already provided as a symbolic `AsynchronousGraph`.
        Alternatively, a "raw" `BooleanNetwork` can be provided and the symbolic
        `AsynchronousGraph` is created automatically.
    succession : list[BooleanSpace]
        A list of sequentially nested trap spaces that specify the target.
    strategy: str
        The searching strategy to use to look for driver nodes. Options are
        'internal' (default) and 'all'.
    max_drivers_per_succession_node: int | None = None,
        The maximum number of drivers that will be tested for a succession
        diagram node. If `None`, then a number of drivers up to the size of the
        succession diagram node's stable motif will be tested.
    forbidden_drivers: set[str] | None
        A set of forbidden drivers that will not be overridden for control. If
        `None`, then all nodes are candidates for control.

    Returns
    -------
    list[ControlOverrides]
        A list of controls. Each control is a list of lists of driver sets,
        represented as state dictionaries. Each list item corresponds to a list
        of drivers for the corresponding trap space in the succession.
    """
    if isinstance(bn, BooleanNetwork):
        bn = AsynchronousGraph(bn)

    control_strategies: list[ControlOverrides] = []
    assume_fixed: BooleanSpace = {}
    for ts in succession:
        control_strategies.append(
            find_drivers(
                bn,
                ts,
                strategy=strategy,
                assume_fixed=assume_fixed,
                max_drivers_per_succession_node=max_drivers_per_succession_node,
                forbidden_drivers=forbidden_drivers,
            )
        )
        ldoi = percolate_space(bn, ts | assume_fixed)
        assume_fixed.update(ldoi)

    return control_strategies


def find_drivers(
    bn: BooleanNetwork | AsynchronousGraph,
    target_trap_space: BooleanSpace,
    strategy: str = "internal",
    assume_fixed: BooleanSpace | None = None,
    max_drivers_per_succession_node: int | None = None,
    forbidden_drivers: set[str] | None = None,
) -> ControlOverrides:
    """
    Finds drives of a given target trap space

    Generally, it is not necessary to call this function directly, as it is
    automatically invoked by
    :func:`drivers_of_succession<biobalm.control.drivers_of_succession>`, which in
    turn is invoked by
    :func:`succession_control<biobalm.control.succession_control>`. It is primarily
    provided in the public API for testing and benchmarking purposes.

    Parameters
    ----------
    bn : BooleanNetwork | AsynchronousGraph
        The network to analyze, which contains the Boolean update functions. Ideally,
        the network should be already provided as a symbolic `AsynchronousGraph`.
        Alternatively, a "raw" `BooleanNetwork` can be provided and the symbolic
        `AsynchronousGraph` is created automatically.
    target_trap_space : BooleanSpace
        The trap space we want to find drivers for.
    strategy: str
        The searching strategy to use to look for driver nodes. Options are
        'internal' (default) and 'all'.
    assume_fixed: dict[str,int] | None
        A dictionary of fixed variables that should be assumed to be fixed.
    max_drivers_per_succession_node: int | None = None,
        The maximum number of drivers that will be tested for a succession
        diagram node. If `None`, then a number of drivers up to the size of the
        succession diagram node's stable motif will be tested.
    forbidden_drivers: set[str] | None
        A set of forbidden drivers that will not be overridden for control. If
        `None`, then all nodes are candidates for control.

    Returns
    -------
    ControlOverrides
        A list of internal driver sets, represented as state dictionaries. If
        empty, then no drivers are found. This can happen if
        `max_drivers_per_succession_node` is not `None`, or if all controls
        require nodes in `forbidden_drivers`.
    """
    if isinstance(bn, BooleanNetwork):
        bn = AsynchronousGraph(bn)

    if assume_fixed is None:
        assume_fixed = {}
    if forbidden_drivers is None:
        forbidden_drivers = set()

    target_trap_space_inner = {
        k: v for k, v in target_trap_space.items() if k not in assume_fixed
    }

    if strategy == "internal":
        driver_pool = set(target_trap_space_inner) - forbidden_drivers
    elif strategy == "all":
        driver_pool = set(bn.network_variable_names()) - forbidden_drivers
    else:
        raise ValueError("Unknown driver search strategy")

    if max_drivers_per_succession_node is None:
        max_drivers_per_succession_node = len(target_trap_space_inner)

    drivers: ControlOverrides = []
    for driver_set_size in range(max_drivers_per_succession_node + 1):
        for driver_set in combinations(driver_pool, driver_set_size):
            if any(set(d) <= set(driver_set) for d in drivers):
                continue

            if strategy == "internal":
                driver_dict: BooleanSpace = {
                    k: cast(Literal[0, 1], target_trap_space_inner[k])
                    for k in driver_set
                }
                ldoi = percolate_space(bn, driver_dict | assume_fixed)
                if target_trap_space.items() <= ldoi.items():
                    drivers.append(driver_dict)
            elif strategy == "all":
                for vals in product([0, 1], repeat=driver_set_size):
                    driver_dict = {
                        driver: cast(Literal[0, 1], value)
                        for driver, value in zip(driver_set, vals)
                    }
                    ldoi = percolate_space(bn, driver_dict | assume_fixed)
                    if target_trap_space.items() <= ldoi.items():
                        drivers.append(driver_dict)
    return drivers


def controls_are_equal(a: ControlOverrides, b: ControlOverrides) -> bool:
    """
    Determine if two :class:`ControlOverrides<biobalm.types.ControlOverrides>`
    objects are equal.

    Two `ControlOverrides` objects are equal if they contain the same
    :class:`BooleanSpace<biobalm.types.BooleanSpace>` objects, regardless of their
    ordering.

    Parameters
    ----------
    a : ControlOverrides
        First :class:`ControlOverrides<biobalm.types.ControlOverrides>` object for
        comparison.
    b : ControlOverrides
        Second :class:`ControlOverrides<biobalm.types.ControlOverrides>` object for
        comparison.

    Returns
    -------
    bool
        Returns `True` if the two
        :class:`ControlOverrides<biobalm.types.ControlOverrides>` objects are
        equal.
    """
    return set(frozenset(x.items()) for x in a) == set(frozenset(x.items()) for x in b)
