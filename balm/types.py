from typing import Literal, TypeAlias, TypedDict

import networkx as nx  # type: ignore

BooleanSpace: TypeAlias = dict[str, Literal[0, 1]]
SubspaceSuccession: TypeAlias = list[BooleanSpace]  # sequence of stable motifs
ControlOverrides: TypeAlias = list[
    BooleanSpace
]  # ways of locking in an individual stable motif


class SuccessionDiagramState(TypedDict):
    network_rules: str
    petri_net: nx.DiGraph
    nfvs: list[str] | None
    dag: nx.DiGraph
    node_indices: dict[int, int]
