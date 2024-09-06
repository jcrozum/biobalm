from __future__ import annotations

from typing import Literal, TypeAlias, TypedDict

import networkx as nx  # type: ignore
import biodivine_aeon as ba

BooleanSpace: TypeAlias = dict[str, Literal[0, 1]]
"""Type alias for `dict[str, Literal[0, 1]]`. Represents a Boolean subspace, which is defined by a set of fixed node values."""
SubspaceSuccession: TypeAlias = list[BooleanSpace]
"""Type alias for `list[BooleanSpace]`. Represents a sequence of subspaces, usually nested."""
ControlOverrides: TypeAlias = list[BooleanSpace]
"""Type alias for `list[BooleanSpace]`. Represents a sequence of subspaces, that percolate to targets."""


class SuccessionDiagramState(TypedDict):
    """
    A `TypedDict` class that stores the state of a succession diagram (see :class:`biobalm.SuccessionDiagram`).
    """

    network_rules: str
    """
    The network rules as an `.aeon` formatted string.
    """

    petri_net: nx.DiGraph
    """
    The Petri net representation of the network rules (see :mod:`biobalm.petri_net_translation`).
    """

    nfvs: list[str] | None
    """
    The negative feedback vertex set used for attractor detection.
    """

    dag: nx.DiGraph
    """
    The directed acyclic graph representation of the succession diagram structure, including
    the :class:`biobalm.types.NodeData` of each node.
    """

    node_indices: dict[int, int]
    """
    A dictionary mapping subspace keys to their positions in the succession
    diagram (see :func:`biobalm.space_utils.space_unique_key`).
    """

    config: SuccessionDiagramConfiguration
    """
    "Global" configuration of a succession diagram.
    """


class NodeData(TypedDict):
    """
    A `TypedDict` class that stores the data of a succession diagram node (see :class:`biobalm.SuccessionDiagram`).

    Returned from :func:`biobalm.SuccessionDiagram.node_data`.
    However, this class is not directly used at runtime, and only exists for static type-checking.
    Instead, at runtime, an untyped dictionary is used because that is what is returned by `networkx.DiGraph.nodes(data=True)`.
    """

    depth: int
    """
    The depth of the node in the succession diagram rooted directed acyclic
    graph. In cases of multiple paths from the root to the reference node,
    the longest is used.
    """

    space: BooleanSpace
    """
    The sub-space that the node represents (this subspace will always be a
    trap space).
    """

    expanded: bool
    """
    Whether the node has been expanded yet or not (`biobalm` builds the
    succession diagram lazily).

    If `expanded=False`, the node *must not* have any successor nodes.
    If `expanded=True`, the node must have all its successors computed
    and included in the succession diagram.
    """

    percolated_network: ba.BooleanNetwork | None
    """
    The AEON `BooleanNetwork` that has variables fixed and percolated
    according to the node `space`. Constant variables are then fully
    eliminated.

    If `None`, this has not been computed yet (`SuccessionDiagram.network` can be
    often used instead).
    """

    percolated_petri_net: nx.DiGraph | None
    """
    The Petri net representation of the network rules percolated to the
    node's sub-space (i.e. a Petri net encoding of the `percolated_network`).
    Constant variables are fully eliminated.

    If `None`, this has not been computed yet (`SuccessionDiagram.petri_net` can be
    often used instead).
    """

    percolated_nfvs: list[str] | None
    """
    An NFVS of the `percolated_network`.

    Note that this contains no variables that are fixed by the node `space`.
    Also, there is no guarantee that this is a subset of `SuccessionDiagram.nfvs`.

    If `None`, this has not been computed yet (`SuccessionDiagram.nfvs` can be
    often used instead).
    """

    attractor_candidates: list[BooleanSpace] | None
    """
    List of states that cover all network attractors in the node's sub-space
    (excluding child nodes, if known).

    Note that each complex/cyclic attractor can be covered by more than
    one state, but each attractor has to be covered by at least one state. Furthermore,
    outside of minimal trap spaces, some candidate states can cover zero attractors.
    That is, even if the candidate set is not empty, it does not guarantee that
    the non-minimal sub-space contains attractors.

    If `None`, these have not been computed, and the node may or may not
    have associated attractors.
    """

    attractor_seeds: list[BooleanSpace] | None
    """
    List of states that one-to-one correspond to network attractors in the node's
    sub-space (excluding child nodes, if known).

    This is very similar to `attractor_candidates`, but here, it is guaranteed
    that each attractor is represented by exactly one state.

    If `None`, these have not been computed, and the node may or may not
    have associated attractors.
    """

    attractor_sets: list[ba.VertexSet] | None
    """
    List of attractors that are present in the node's sub-space (excluding
    child nodes, if known).

    Each attractor is represented symbolically using `biodivine_aeon.VertexSet`.

    Note that for fixed-point attractors, the attractor set is effectively
    equivalent to the attractor seed. However, for complex attractors, this
    set containes *all* the attractor states as opposed to just one. Hence it
    is harder to compute, but can be used to analyze things like the average
    value of each variable in the attractor states, or the presence of
    particular oscillation patterns.

    If `None`, these have not been computed, and the node may or may not
    have associated attractors.
    """

    parent_node: int | None
    """
    The ID of one of the predecessor nodes through which this node was first
    discovered, or `None` if the node was discovered through other means
    (usually either the root, or a minimal trap space that was added without
    expanding its predecessors).
    """

    skipped: bool | None
    """
    If `True`, indicates that the successors of this node are not the normal
    maximal trap spaces, but rather some smaller trap spaces deeper in the
    succession diagram. The only requirement is that all minimal trap spaces
    that are reachable from this "skip node" are still reachable through
    this modified set of successors (in the most basic form, the successors
    can be the minimal trap spaces themselves).

    For skip nodes, attractor detection still works, but can over-count
    motif-avoidant attractors assuming the attractor intersects multiple
    skip nodes. In other words, each attractor is still found, but it is
    not necessarily true that its smallest trap space is discovered.
    Nevertheless, if the network has no motif-avoidant attractors, the
    result is always the same regardless of how many skip nodes are used.

    A skip node is considered "expanded", since its outgoing edges are
    computed.
    """


class SuccessionDiagramConfiguration(TypedDict):
    """
    Describes the configuration options of a `SuccessionDiagram`.

    Use :meth:`SuccessionDiagram.default_config` to create a
    configuration dictionary pre-populated with default values.
    """

    debug: bool
    """
    If `True`, the `SuccessionDiagram` will print messages
    describing the progress of the running operations.

    [Default: False]
    """

    max_motifs_per_node: int
    """
    Limit on the number of stable motifs explored for one node of a succession
    diagram. If this limit is exceeded during node expansion, a `RuntimeError`
    is raised and the node remains unexpanded.

    This limit is in place mainly to avoid surprising out of memory errors,
    because currently there is no logging mechanism that would report the
    number of stable motifs gradually.

    [Default: 100_000]
    """

    nfvs_size_threshold: int
    """
    For networks larger than this threshold, we only run FVS detection
    instead of NFVS detection. This is still correct, but can produce
    a larger node set.


    There is a trade-off between the speed gains from a smaller node set
    to consider and the cost of determining which FVS nodes only
    intersect negative cycles to find an NFVS subset. Typically,
    for smaller networks, the trade-off is in favor of
    computing a smaller NFVS.
    """

    pint_goal_size_limit: int
    """
    Pint is called using command line and can only accept a limited number of arguments.
    This limit is applied when constructing goals to avoid reaching this argument limit.

    The limit currently applies to the total number of literals that can be used to
    represent a goal.

    The default value was empirically tested as safe on Debian linux, but other operating
    systems may need a different limit to stay safe. Nevertheless, this should not be
    an issue on smaller/simpler networks.
    """

    attractor_candidates_limit: int
    """
    If more than `attractor_candidates_limit` states are produced during the
    attractor detection process, then the process fails with a `RuntimeError`.
    This is mainly to avoid out-of-memory errors or crashing `clingo`.
    """

    retained_set_optimization_threshold: int
    """
    If there are more than this amount of attractor candidates, the attractor
    detection process will try to optimize the retained set using ASP (if enabled).
    """

    minimum_simulation_budget: int
    """
    The minimum number of simulation steps per network variable that is guaranteed to be
    spent on eliminating attractor candidate states. That is, if budget is `1_000` and
    the network has `200` variables, the total number of allowed
    simulation steps is `200_000`.

    Note that this is a budget that applies to all candidates collectively. So if the number
    of candidates is larger, the number of steps per candidate is proportionally smaller.
    However, this budget only applies when simulation has not been able to make progress
    in the recent round. That is, if simulation has actively eliminated some candidates in
    the recent round, it will still continue regardless of the budget limit.
    """
