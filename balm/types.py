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
    """A `TypedDict` class that stores the state of a succession diagram (see :class:`balm.SuccessionDiagram`).

    Parameters
    ----------
    network_rules : str
        The network rules as a string.
    petri_net : networkx.DiGraph
        The Petri net representation of the network rules.
    nfvs : list[str] | None
        The negative feedback vertex set used for attractor detection.
    dag : networkx.DiGraph
        The directed acyclic graph (DAG) representation of the succession
        diagram structure.
    node_indices : dict[int, int]
        A dictionary mapping subspace keys to their positions in the succession
        diagram (see :func:`balm.space_utils.space_unique_key`).
    """

    network_rules: str
    petri_net: nx.DiGraph
    nfvs: list[str] | None
    dag: nx.DiGraph
    node_indices: dict[int, int]
