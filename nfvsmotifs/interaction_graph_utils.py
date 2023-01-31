from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from biodivine_aeon import VariableId # type: ignore
    from typing import Union, Optional

from biodivine_aeon import BooleanNetwork, RegulatoryGraph
from networkx import DiGraph # type: ignore

from typing import List, Set
from networkx.algorithms import bipartite # type: ignore

from pyeda.boolalg.expr import expr # type: ignore
from pyeda.boolalg.bdd import expr2bdd, bddvar # type: ignore

from nfvsmotifs.SignedGraph import SignedGraph

"""
A python package for approximating minimum feedback vertex sets
"""
from nfvsmotifs.FVSpython3 import FVS as FVS

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
    network = network.infer_regulatory_graph()

    # Convert AEON `RegulatoryGraph` to a generic `DiGraph` that we can use later.
    rg = network.graph()
    ig = DiGraph()

    for var in rg.variables():
        ig.add_node(rg.get_variable_name(var))

    for reg in rg.regulations():
        if not reg['observable']:
            # In general, a fully specified network should only contain 
            # observable (i.e. essential) regulations.
            raise Exception("Unreachable: You are using this on partially specified networks, aren't you?")
        source = rg.get_variable_name(reg['source'])
        target = rg.get_variable_name(reg['target'])
        sign = "?"        
        if 'monotonicity' in reg:
            # Monotonic regulation---add only one edge.
            if reg['monotonicity'] == 'activation':
                sign = "+"
            elif reg['monotonicity'] == 'inhibition':
                sign = "-"
            else:
                raise Exception(f"Unreachable: unknown monotonicity {reg['monotonicity']}.")
        ig.add_edge(source, target, sign=sign)
    return ig

def _digraph_to_regulatory_graph(graph: DiGraph) -> RegulatoryGraph:
    """
        A helper method to transform between a "signed digraph" and AEON's `RegulatoryGraph`.

        All nodes of the digraph should have string identifiers. Edges can be optionally 
        annotated with a `sign` value `"+"`, `"-"` or `"?"` (however, `"?"` is treated the 
        same as a missing `sign` annotation).
    """
    rg = RegulatoryGraph(list(graph.nodes()))
    for edge in graph.edges():
        edge_data = graph.get_edge_data(edge[0], edge[1])
        monotonicity = None
        if 'sign' in edge_data:
            sign = edge_data['sign']
            if sign == '+':
                monotonicity = 'activation'
            elif sign == '-':
                monotonicity = 'inhibition'
            elif sign != '?':
                raise Exception(f"Unknown monotonicity sign: '{sign}'. Expected '+'/'-'/'?'")
        rg.add_regulation({
            'source': edge[0],
            'target': edge[1],
            'observable': True, # For now, observability is not in the graph.
            'monotonicity': monotonicity
        })

    return rg


def feedback_vertex_set(
    network: Union[BooleanNetwork, RegulatoryGraph, DiGraph], 
    parity: Optional[str] = None, 
    subgraph: Optional[list[Union[str, VariableId]]] = None
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
    if type(network) == BooleanNetwork:
        network = network.graph()
    if type(network) == DiGraph:
        network = _digraph_to_regulatory_graph(network)
    fvs = network.feedback_vertex_set(parity=parity, restriction=subgraph)
    return [network.get_variable_name(x) for x in fvs]

def independent_cycles(
    network: Union[BooleanNetwork, RegulatoryGraph],
    parity: Optional[str] = None,
    subgraph: Optional[list[Union[str, VariableId]]] = None
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
    if type(network) == BooleanNetwork:
        network = network.graph()
    if type(network) == DiGraph:
        network = _digraph_to_regulatory_graph(network)
    ic = network.independent_cycles(parity=parity, restriction=subgraph)
    return [[network.get_variable_name(x) for x in cycle] for cycle in ic]


def find_minimum_NFVS(network: BooleanNetwork) -> list[str]:

    """
    BDD variables
    """
    bdd_vars = {}

    """
    BDDs of Boolean functions
    """
    bdd_funs = {}

    """
    Negative feedback vertex set
    """
    U_neg = []

    """
    List of source nodes
    """

    source_nodes = []

    """
    Node to input nodes
    """

    INx = {}

    nodes = []

    for variable in network.variables():
        var_name = network.get_variable_name(variable)
        function = network.get_update_function(variable)

        nodes.append(var_name)

        if function.strip() == var_name:
            source_nodes.append(var_name)

        fx = expr(function.replace("!", "~"))

        INx[var_name] = fx.support # list of nodes appearing in Boolean function fx

        vx = bddvar(var_name)
        fx = expr2bdd(expr(fx))

        bdd_vars[var_name] = vx
        bdd_funs[var_name] = fx

            
    """Build the unsigned and signed interaction graphs"""
    u_ig = DiGraph()
    s_ig = SignedGraph(nodes)
    
    for x in nodes:
        u_ig.add_node(x)

        fx = bdd_funs[x]

        for y in INx[x]:
            is_actual_arc = False

            vy = bdd_vars[str(y)]
            
            fx_res_vy_0 = fx.restrict({vy: 0})
            fx_res_vy_1 = fx.restrict({vy: 1})

            pos_arc = ~fx_res_vy_0 & fx_res_vy_1
            neg_arc = fx_res_vy_0 & ~fx_res_vy_1

            if pos_arc.is_one() or pos_arc.satisfy_one():
                # a positive arc with weight = 1
                s_ig.set_edge(str(y), str(x), 1)
                is_actual_arc = True

            if neg_arc.is_one() or neg_arc.satisfy_one():
                # a negative arc with weight = -1
                s_ig.set_edge(str(y), str(x), -1)
                is_actual_arc = True

            if is_actual_arc == True:
                u_ig.add_edge(str(y), str(x))


    """First, find feedback vertex set"""
    U = FVS.FVS(u_ig)

    U = list(set(U) - set(source_nodes))

    """Second, filter feedback vertex set to get an negative feedback vertex set"""
    U_neg = s_ig.get_self_negative_loops()
    U_candidate = []

    for v in U:
        if not v in U_neg:
            U_candidate.append(v)

    for v in U_neg:
        s_ig.remove_vertex(v)

    for v in source_nodes:
        s_ig.remove_vertex(v)
    
    while not is_no_negative_cycle(s_ig):
        v = select_by_negative_degree(s_ig, U_candidate)

        if len(v) == 0:
            break
    
        U_candidate.remove(v)

        U_neg.append(v)
        
        s_ig.remove_vertex(v)

    
    return U_neg


def is_no_negative_cycle(s_ig: SignedGraph) -> bool:
    udGraph = s_ig.convert_to_undirected_graph()
    return bipartite.is_bipartite(udGraph)


def select_by_negative_degree(s_ig: SignedGraph, U_candidate: list[str]) -> str:
    v_selected = ""
    max_neg_deg = -1

    for v in U_candidate:
        v_neg_deg = s_ig.get_negative_degree(v)

        if v_neg_deg >= max_neg_deg:
            v_selected = v
            max_neg_deg = v_neg_deg

    return v_selected