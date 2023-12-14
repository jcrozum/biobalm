from typing import Literal, TypeAlias

BooleanSpace: TypeAlias = dict[str, Literal[0, 1]]
SubspaceSuccession: TypeAlias = list[BooleanSpace]  # sequence of stable motifs
ControlOverrides: TypeAlias = list[
    BooleanSpace
]  # ways of locking in an individual stable motif
