from typing import Literal

space_type = dict[str, Literal[0, 1]]
SuccessionType = list[space_type]  # sequence of stable motifs
ControlType = list[space_type]  # ways of locking in an individual stable motif
