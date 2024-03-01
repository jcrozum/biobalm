from __future__ import annotations

from typing import Literal, TypeAlias, TypedDict

import networkx as nx  # type: ignore

BooleanSpace: TypeAlias = dict[str, Literal[0, 1]]
"""Type alias for `dict[str, Literal[0, 1]]`. Represents a Boolean subspace, which is defined by a set of fixed node values."""
SubspaceSuccession: TypeAlias = list[BooleanSpace]
"""Type alias for `list[BooleanSpace]`. Represents a sequence of subspaces, usually nested."""
ControlOverrides: TypeAlias = list[BooleanSpace]
"""Type alias for `list[BooleanSpace]`. Represents a sequence of subspaces, that percolate to targets."""


class SuccessionDiagramState(TypedDict):
    """
    A `TypedDict` class that stores the state of a succession diagram (see :class:`balm.SuccessionDiagram`).
    """

    network_rules: str
    """
    The network rules as a string.
    """

    petri_net: nx.DiGraph
    """
    The Petri net representation of the network rules.
    """

    nfvs: list[str] | None
    """
    The negative feedback vertex set used for attractor detection.
    """

    dag: nx.DiGraph
    """
    The directed acyclic graph representation of the succession diagram structure.
    """

    node_indices: dict[int, int]
    """
    A dictionary mapping subspace keys to their positions in the succession
    diagram (see :func:`balm.space_utils.space_unique_key`).
    """


class NodeData(TypedDict):
    """
    A `TypedDict` class that stores the data of a succession diagram node (see :class:`balm.SuccessionDiagram`).

    Returned from :func:`balm.SuccessionDiagram.node_data`.
    However, this class is not directly used at runtime, and only exists for static type-checking.
    Instead, at runtime, an untyped dictionary is used because that is what is returned by `networkx.DiGraph.nodes(data=True)`.
    """

    depth: int
    """
    The depth of the node in the succession diagram rooted directed acyclic
    graph. In cases of multiple paths from the root to the reference node,
    the longest is used.
    """

    attractors: list[BooleanSpace] | None
    """
    Attractor seed states for the node. If `None`, these have not been
    computed, and the node may or may not have associated attractors.
    """

    petri_net: nx.DiGraph | None
    """
    The Petri net representation of the node. If `None`, this have not been
    computed yet.
    """
    space: BooleanSpace
    """
    The sub-space that the node represents (this subspace will always be a
    trap space).
    """

    expanded: bool
    """
    Whether the node has been expanded yet or not (`balm` builds the
    succession diagram lazily).
    """
