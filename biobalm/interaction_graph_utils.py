"""
A utility module which provides several methods for manipulating the influence graph
of a Boolean network.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

if TYPE_CHECKING:
    from typing import Any, Sequence

    from biodivine_aeon import Regulation, VariableId

from typing import cast

from biodivine_aeon import BooleanNetwork, RegulatoryGraph, SignType, SymbolicContext
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
        edges are considered. Edges with no monotonicity are counted as both
        positive and negative. Otherwise, all cycles are considered.
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

    Examples
    --------
    >>> import biobalm
    >>> from biobalm.interaction_graph_utils import feedback_vertex_set
    >>> sd = biobalm.SuccessionDiagram.from_rules(\"\"\"
    ...     A, B
    ...     B, A
    ...     C, D
    ...     D, !C\"\"\")
    >>> feedback_vertex_set(sd.network)
    ['A', 'C']
    >>> feedback_vertex_set(sd.network, parity="positive")
    ['A']
    >>> feedback_vertex_set(sd.network, parity="negative")
    ['C']
    """
    if isinstance(network, DiGraph):
        network = _digraph_to_regulatory_graph(network)
    assert isinstance(network, RegulatoryGraph)
    fvs = network.feedback_vertex_set(parity, subgraph)
    return sorted([network.get_variable_name(x) for x in fvs])


def cleanup_network(network: BooleanNetwork) -> BooleanNetwork:
    """
    Prepare a `BooleanNetwork` object for use in a `SuccessionDiagram`. This
    mainly ensures that the network does not use parameters and fixes all static
    constraints to ensure that they are actually correct.

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

    return network.infer_valid_graph()


def source_SCCs(bn: BooleanNetwork) -> list[list[str]]:
    """
    Find source SCCs of the given `BooleanNetwork`.

    Here, SCC stands for "strongly connected component". An SCC is a source SCC
    if it has no incoming edges.

    Parameters
    ----------
    bn : BooleanNetwork
        The Boolean network to be examined.

    Returns
    -------
    list[list[str]]
        The list of source SCCs.
    """
    result: list[list[str]] = []
    for scc in bn.strongly_connected_components():
        scc_list = sorted(scc)
        if bn.backward_reachable(scc_list) == scc:
            scc_names = [bn.get_variable_name(var) for var in scc_list]
            result.append(scc_names)

    return sorted(result)


def source_nodes(
    network: BooleanNetwork, ctx: SymbolicContext | None = None
) -> list[str]:
    """
    Identify the source nodes of a given `BooleanNetwork`.

    Return the source nodes of a `BooleanNetwork`. That is, variables whose
    value cannot change, but is not fixed to a `true`/`false` constant.

    Note that this internally uses BDD translation to detect identity functions
    semantically rather than syntactically. If you already have a
    `SymbolicContext` for the given `network` available, you can supply it as
    the second argument.

    Parameters
    ----------
    network : BooleanNetwork
        The Boolean network to be examined.
    ctx : SymbolicContext
        The context used to translate the network to BDDs. A
        `biodivine_aeon.SymbolicContext` object.

    Returns
    -------
    list[str]
        The list of source nodes.
    """
    if ctx is None:
        ctx = SymbolicContext(network)

    result: list[str] = []
    for var in network.variable_names():
        update_function = network.get_update_function(var)
        if update_function is None:
            # This is an input variable with unspecified update
            # (this defaults to identity in most tools).
            assert len(network.predecessors(var)) == 0
            result.append(var)
        else:
            fn_bdd = ctx.mk_update_function(update_function)
            var_bdd = ctx.mk_network_variable(var)
            if fn_bdd == var_bdd:
                result.append(network.get_variable_name(var))

    return result
