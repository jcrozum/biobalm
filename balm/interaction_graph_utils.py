from __future__ import annotations

from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from biodivine_aeon import VariableId, Regulation
    from typing import Any

from typing import cast
from biodivine_aeon import BooleanNetwork, RegulatoryGraph, SignType
from networkx import DiGraph  # type: ignore


def infer_signed_interaction_graph(network: BooleanNetwork) -> DiGraph:
    """
    Takes an arbitrary `BooleanNetwork` and extracts a signed interacion graph
    based on the *actual* dependencies of the network's update functions.

    In particular, this eliminates things like unused variables and it will also
    infer correct monotonicity of inputs. Importantly, this is performed using
    BDDs, so it is not merely a syntactic transformation (but general BDD
    limitations apply).

    In the resulting digraph, we use edge attribute `sign="+"` or `sign="-"`
    to differentiate between positive and negative edges. An iteraction that
    is not monotonic is represented as `sign="?"`.
    """
    # Since the model files can be incorrect, we have to infer the regulations
    # from the network functions directly. Note that this step uses BDDs, so in the
    # unlikely event that the update functions are too complex to represent as BDDs,
    # this will crash on an out-of-memory error, or just run for a very long time.
    network = network.infer_valid_graph()

    # Convert AEON `RegulatoryGraph` to a generic `DiGraph` that we can use later.
    rg = network.to_graph()
    ig = DiGraph()

    for var in rg.variables():
        ig.add_node(rg.get_variable_name(var))  # type: ignore

    for reg in rg.regulations():
        if not reg["essential"]:
            # In general, a fully specified network should only contain
            # essential regulations after `infer_valid_graph`.
            raise Exception(
                "Unreachable: You are using this on partially specified networks, aren't you?"
            )
        source = rg.get_variable_name(reg["source"])
        target = rg.get_variable_name(reg["target"])
        sign = reg["sign"] if reg["sign"] is not None else "?"
        ig.add_edge(source, target, sign=sign)  # type: ignore
    return ig


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
        reg: Regulation[str] = {
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
    network: BooleanNetwork | RegulatoryGraph | DiGraph,
    parity: Literal["positive", "negative"] | None = None,
    subgraph: list[str] | list[VariableId] | None = None,
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
    if isinstance(network, BooleanNetwork):
        network = network.to_graph()
    if isinstance(network, DiGraph):
        network = _digraph_to_regulatory_graph(network)
    assert isinstance(network, RegulatoryGraph)
    fvs = network.feedback_vertex_set(parity=parity, subgraph=subgraph)
    return sorted([network.get_variable_name(x) for x in fvs])


def independent_cycles(
    network: BooleanNetwork | RegulatoryGraph,
    parity: Literal["positive", "negative"] | None = None,
    subgraph: list[str] | list[VariableId] | None = None,
) -> list[list[str]]:
    """
    Compute an approximately maximal set of independent cycles of
    a `BooleanNetwork`, `RegulatoryGraph` or a `DiGraph` with optional `sign` annotations
    on its edges. The result is guaranteed to be a set of independent cycles, but it may
    not be maximal.

    There are two optional parameters:

     - `parity`: Can be either `positive` or `negative`. If parity is specified, only cycles of
        the specified parity are considered (e.g. if `pairty='negative'`, there can still be positive
        cycles in the graph which are not covered by this independent cycle set).
     - `subgraph`: A list of network variables (either string names or AEON `VariableId` objects
        are fine). If given, the result is restricted to the sub-graph induced by these network
        nodes.

    The result is a list of cycles, such that each cycle is a list of variable names in the order
    in which they appear on the cycle. The cycles are sorted by increasing length.

    In general, the method should be deterministic (the same pseudo-optimal cycles are returned
    every time). However, while I believe the sorting should be stable too, please treat the
    order of returned cycles with caution :)
    """
    if isinstance(network, BooleanNetwork):
        network = network.to_graph()

    # this should never happen, but it's easy enough to convert
    if isinstance(network, DiGraph):
        network = _digraph_to_regulatory_graph(network)

    ic = network.independent_cycles(parity=parity, subgraph=subgraph)
    return [[network.get_variable_name(x) for x in cycle] for cycle in ic]
