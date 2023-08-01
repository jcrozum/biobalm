"""
This is an internal module which contains the implementations of the more
involved algorithms related to `SuccessionDiagram` construction and
manipulation. The idea is that any implementation that requires some helper
methods or data structures should be present here instead of
`SuccessionDiagram.py` to avoid mixing of code used by different algorithms and
bloating the core succession diagram implementation.

The algorithms within this module here should only use the public API of the
`SuccessionDiagram` to avoid violating its invariants. The actual algorithms are
then "re-exported" as a public methods of the `SuccessionDiagram` class.
"""
