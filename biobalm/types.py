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
