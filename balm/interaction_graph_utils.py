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
    A helper method to convert a `networkx.DiGraph` to a `biodivine_aeon.RegulatoryGraph`.

    All nodes of the digraph should have string identifiers. Edges can be optionally
    annotated with a `sign` value `"+"`, `"-"` or `"?"` (however, `"?"` is treated the
    same as a missing `sign` annotation).

    Parameters
    ----------
    graph : DiGraph
        A `networkx.DiGraph` with string identifiers as nodes and optional
        `sign` annotations for edges.

    Returns
    -------
    RegulatoryGraph
        A `biodivine_aeon.RegulatoryGraph` with the same nodes and edges as `graph`.
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
    Compute an approximately minimal feedback vertex set (FVS).

    Uses the implementation of `biodivine_aeon`.

    An FVS is a set of nodes in a network whose removal from the network
    results in an acyclic network. An FVS is minimal if it does not contain any
    FVS as a proper subset. This method uses a heuristic approach to attempt to find a
    minimal FVS. The result is guaranteed to be a feedback vertex set, but it might not be minimal.

    A paritiy can be specified. If the parity is specified, only cycles of the
    specified parity are considered (e.g. if `parity='negative'`, there can still be positive
    cycles in the graph not covered by the returned feedback vertex set).

    The method is deterministic (the same pseudo-optimal FVS is returned every time).

    Parameters
    ----------
    network : RegulatoryGraph | DiGraph
        Network to consider. If a `DiGraph` is given, it is converted to a `RegulatoryGraph`.
    parity : Literal["positive";, "negative"] | None, optional
        Must be `"positive"`, `"negative"`, or `None` (default). If
        `"positive"`, only cycles with an even number of negative edges are
        considered. If `"negative"`, only cycles with an odd number of negative
        edges are considered. Otherwise, all cycles are considered.
    subgraph : Sequence[str  |  VariableId] | None, optional
        A list of network variables (either string names or AEON `VariableId`
        objects are fine). If given, the result is the FVS of the sub-graph induced by
        these network nodes.

    Returns
    -------
    list[str]
        A list of node names representing the smallest FVS found during the
        search. Sorted in the same order as in the input network (usually
        lexicographically).
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

    Parameters
    ----------
    network : BooleanNetwork
        The Boolean network to be prepared.

    Returns
    -------
    BooleanNetwork
        The prepared network.
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
