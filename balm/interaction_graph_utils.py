from __future__ import annotations

from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from typing import Any, Sequence

    from biodivine_aeon import Regulation, VariableId

import copy
from typing import cast

from biodivine_aeon import BooleanNetwork, RegulatoryGraph, SignType
from networkx import DiGraph  # type: ignore


def _digraph_to_regulatory_graph(graph: DiGraph) -> RegulatoryGraph:
    """
    A helper method to transform between a "signed digraph" and AEON's `RegulatoryGraph`.

    All nodes of the digraph should have string identifiers. Edges can be optionally
    annotated with a `sign` value `"+"`, `"-"` or `"?"` (however, `"?"` is treated the
    same as a missing `sign` annotation).
    """
    rg = RegulatoryGraph(list(graph.nodes()))  # type: ignore
    for edge in graph.edges():  # type: ignore
        edge_data: dict[Any, Any] = graph.get_edge_data(edge[0], edge[1])  # type: ignore
        reg: Regulation = {
            "source": edge[0],
            "target": edge[1],
            "essential": True,
            "sign": None,
        }
        if "sign" in edge_data and edge_data["sign"] in ["+", "-"]:
            reg["sign"] = cast(SignType, edge_data["sign"])

        rg.add_regulation(reg)

    return rg


def feedback_vertex_set(
    network: RegulatoryGraph | DiGraph,
    parity: Literal["positive", "negative"] | None = None,
    subgraph: Sequence[str | VariableId] | None = None,
) -> list[str]:
    """
    Compute an approximately minimal feedback vertex set (FVS) of
    a `BooleanNetwork`, `RegulatoryGraph` or a `DiGraph` with optional `sign` annotations
    on its edges. The result is guaranteed to be a feedback vertex set, but it might not be minimal.

    There are two optional parameters:

     - `parity`: Can be either `positive` or 'negative'. If parity is specified, only cycles of the
        specified parity are considered (e.g. if `parity='negative'`, there can still be positive
        cycles in the graph not covered by the returned feedback vertex set).
     - `subgraph`: A list of network variables (either string names or AEON `VariableId`
        objects are fine). If given, the result is the FVS of the sub-graph induced by
        these network nodes.

    The result is a list of variable names representing the FVS. The variables are always sorted
    based on the order in which they appear in the network (which is typically lexicographic,
    if the network is loaded from a file and not made "by hand").

    The method should be deterministic (the same pseudo-optimal FVS is returned every time).
    """
    if isinstance(network, DiGraph):
        network = _digraph_to_regulatory_graph(network)
    assert isinstance(network, RegulatoryGraph)
    fvs = network.feedback_vertex_set(parity, subgraph)
    return sorted([network.get_variable_name(x) for x in fvs])


def cleanup_network(network: BooleanNetwork) -> BooleanNetwork:
    """
    Prepare a `BooleanNetwork` object for use in a `SuccessionDiagram`. This mainly
    checks that the network has no parameters and removes any constraints that could
    add additional overhead to symbolic manipulation.
    """

    assert (
        network.explicit_parameter_count() == 0
    ), f"Parametrized networks are not supported. Found parameters: {network.explicit_parameter_names()}."

    # Implicit parameters with no regulators are allowed, since they just reprtesent free inputs
    # and are explicitly handled by the succession diagram.
    non_input_implicit = [
        v for v in network.implicit_parameters() if len(network.predecessors(v)) > 0
    ]
    if len(non_input_implicit) > 0:
        names = [network.get_variable_name(x) for x in non_input_implicit]
        raise AssertionError(
            f"Parametrized networks are not supported. Found implicit parameters: {names}."
        )

    network = copy.copy(network)
    for reg in network.regulations():
        reg["essential"] = False
        reg["sign"] = None
        assert network.ensure_regulation(reg) is not None

    return network
