import os
from biodivine_aeon import BooleanNetwork # type:ignore
from networkx import DiGraph # type:ignore
from nfvsmotifs.interaction_graph_utils import infer_signed_interaction_graph, feedback_vertex_set, independent_cycles, find_minimum_NFVS

def test_ig_inference():
    bn = BooleanNetwork.from_bnet("""
        # Just a normal function.
        b, a | !b
        # Contradiciton on `a` - the regulation should not appear in the result
        # Also, non-monotonic dependence on b and c.
        a, (a & !a) | (b <=> c)
        c, c
    """)
    ig = infer_signed_interaction_graph(bn)

    edges = { edge:ig.get_edge_data(edge[0], edge[1])['sign'] for edge in ig.edges }
    assert len(edges) == 5
    assert edges[('a', 'b')] == "+"
    assert edges[('b', 'b')] == "-"
    assert edges[('b', 'a')] == "?"
    assert edges[('c', 'a')] == "?"
    assert edges[('c', 'c')] == "+"
    assert ('a', 'a') not in edges

# There should be a negative cycle between b_1 and b_2, 
# a positive cycle between d_1 and d_2, and a negative cycle
# between d_1, d_2, and d_3. Other nodes are not on cycles
# except for e, which has a positive self-loop.
CYCLES_BN = BooleanNetwork.from_aeon("""
            a -> c
            b_1 -> b_2
            b_2 -| b_1
            b_2 -> c
            c -> d_2
            c -> e
            d_1 -> d_3
            d_3 -| d_2
            d_2 -> d_1
            d_1 -> d_2
            e -> e
    """)

CYCLES_DIGRAPH = DiGraph()
CYCLES_DIGRAPH.add_nodes_from(["a", "b_1", "b_2", "c", "d_1", "d_2", "d_3", "e"])
CYCLES_DIGRAPH.add_edges_from([
    ("a", "c", {'sign': '+'}),
    ("b_1", "b_2", {'sign': '+'}),
    ("b_2", "b_1", {'sign': '-'}),
    ("b_2", "c", {'sign': '+'}),
    ("c", "d_2", {'sign': '+'}),
    ("c", "e", {'sign': '+'}),
    ("d_1", "d_3", {'sign': '+'}),
    ("d_3", "d_2", {'sign': '-'}),
    ("d_2", "d_1", {'sign': '+'}),
    ("d_1", "d_2", {'sign': '+'}),
    ("e", "e", {'sign': '+'}),
])

def test_fvs():
    fvs = feedback_vertex_set(CYCLES_BN)
    nfvs = feedback_vertex_set(CYCLES_BN, parity='negative')
    pfvs = feedback_vertex_set(CYCLES_BN, parity='positive')

    assert len(fvs) == 3
    assert len(nfvs) == 2
    assert len(pfvs) == 2

    # All of these assertions assume that the greedy algorithm
    # picked a truly minimal FVS, which seems to be the case
    # for a simple graph like this one.

    # "a" and "c" are not on any cycle, so should not be in any FVS
    assert ("a" not in fvs) and ("a" not in nfvs) and ("a" not in pfvs)
    assert ("c" not in fvs) and ("c" not in nfvs) and ("c" not in pfvs)

    # "e" has a positive self-loop, so should be in both fvs and pfvs
    assert ("e" in fvs) and ("e" not in nfvs) and ("e" in pfvs)

    # "b_1" or "b_2" should appear in both fvs and nfvs, but not pfvs
    assert ("b_1" in fvs) or ("b_2" in fvs)
    assert ("b_1" in nfvs) or ("b_2" in nfvs)
    assert ("b_1" not in pfvs) and ("b_2" not in pfvs)
    
    # With "d_*", its a bit more complicated:
    # "d_1" or "d_2" must be in fvs and also pfvs, but in nfvs, "d_3"
    # is also sufficient as the "d_1 --- d_2" cycle is positive.
    assert ("d_1" in fvs) or ("d_2" in fvs)
    assert ("d_1" in nfvs) or ("d_2" in nfvs) or ("d_3" in nfvs)
    assert ("d_1" in pfvs) or ("d_2" in pfvs)
    
    # Check that the `DiGraph` results are the same as `BooleanNetwork` results.
    dg_fvs = feedback_vertex_set(CYCLES_DIGRAPH)
    dg_nfvs = feedback_vertex_set(CYCLES_DIGRAPH, parity='negative')
    dg_pfvs = feedback_vertex_set(CYCLES_DIGRAPH, parity='positive')

    assert fvs == dg_fvs
    assert nfvs == dg_nfvs
    assert pfvs == dg_pfvs

def test_subgraph_fvs():
    # We only keep the two cycles consisting of "d_*". The "b_*" cycle 
    # and "e" self-loop are not considered.
    fvs = feedback_vertex_set(CYCLES_BN, subgraph=["a", "b_1", "d_1", "d_2", "d_3"])
    pfvs = feedback_vertex_set(CYCLES_BN, parity='positive', subgraph=["a", "b_1", "d_1", "d_2", "d_3"])
    nfvs = feedback_vertex_set(CYCLES_BN, parity='negative', subgraph=["a", "b_1", "d_1", "d_2", "d_3"])

    assert len(fvs) == 1
    assert len(pfvs) == 1
    assert len(nfvs) == 1
    assert ("d_1" in fvs) or ("d_2" in fvs)
    assert ("d_1" in nfvs) or ("d_2" in nfvs) or ("d_3" in nfvs)
    assert ("d_1" in pfvs) or ("d_2" in pfvs)

def test_ic():
    ic = independent_cycles(CYCLES_BN)
    n_ic = independent_cycles(CYCLES_BN, parity='negative')
    p_ic = independent_cycles(CYCLES_BN, parity='positive')

    assert len(ic) == 3
    assert len(n_ic) == 2
    assert len(p_ic) == 2
            
    # e is the shortes positive (and overall) cycle, so should be first
    assert ic[0] == ["e"]
    assert p_ic[0] == ["e"]

    # "b_*" is the smallest negative cycle
    assert set(n_ic[0]) == set(["b_1", "b_2"])

    # The second positive cycle is in the shorter "d_*" cycle.
    assert set(p_ic[1]) == set(["d_1", "d_2"])
    # And the second negative cycle is the longer "d_*" cycle.
    assert set(n_ic[1]) == set(["d_1", "d_2", "d_3"])

    # For the general case, both cycles of length two are included.
    # But their order is not guaranteed.
    assert set(ic[1]) == set(["b_1", "b_2"]) or set(ic[2]) == set(["b_1", "b_2"])
    assert set(ic[1]) == set(["d_1", "d_2"]) or set(ic[2]) == set(["d_1", "d_2"])

    # Check that the `DiGraph` results are the same as `BooleanNetwork` results.
    # Note that these are not necessarily entirely equivalent, as the DiGraph 
    # seems to store the nodes/edges in a hashmap, resulting in 
    # not-quite-deterministic ordering and possibly different results (I think?).
    dg_ic = independent_cycles(CYCLES_DIGRAPH)
    dg_n_ic = independent_cycles(CYCLES_DIGRAPH, parity='negative')
    dg_p_ic = independent_cycles(CYCLES_DIGRAPH, parity='positive')

    print(ic)
    print(dg_ic)

    assert [set(x) for x in ic] == [set(x) for x in dg_ic]
    assert [set(x) for x in n_ic] == [set(x) for x in dg_n_ic]
    assert [set(x) for x in p_ic] == [set(x) for x in dg_p_ic]


def test_subgraph_ic():
    # We only keep the two cycles consisting of "d_*". The "b_*" cycle 
    # and "e" self-loop are not considered.
    ic = independent_cycles(CYCLES_BN, subgraph=["a", "b_1", "d_1", "d_2", "d_3"])
    p_ic = independent_cycles(CYCLES_BN, parity='positive', subgraph=["a", "b_1", "d_1", "d_2", "d_3"])
    n_ic = independent_cycles(CYCLES_BN, parity='negative', subgraph=["a", "b_1", "d_1", "d_2", "d_3"])

    assert len(ic) == 1
    assert len(p_ic) == 1
    assert len(n_ic) == 1
    assert set(ic[0]) == set(["d_1", "d_2"]) or set(ic[0]) == set(["d_1", "d_2", "d_3"])
    assert set(p_ic[0]) == set(["d_1", "d_2"])
    assert set(n_ic[0]) == set(["d_1", "d_2", "d_3"])


def test_fvs_accuracy_CASCADE3():
    """
    Compare results of AEON and mtsNFVS on computing an negative feedback vertex set of the CASCADE3 model <https://doi.org/10.3389/fmolb.2020.502573>.
    Note that the result of mtsNFVS is not deterministic.
    """
    file_path = os.getcwd() + "/CASCADE3.bnet"
    bn_real = BooleanNetwork.from_file(file_path)

    nfvs_aeon = feedback_vertex_set(bn_real, parity='negative')
    nfvs_mtsNFVS = find_minimum_NFVS(bn_real)

    assert len(nfvs_mtsNFVS) <= 19 # the result of mtsNFVS is 19
    assert len(nfvs_aeon) <= 22 # the result of AEON is 23
    

def test_fvs_accuracy_SIPC():
    """
    Compare results of AEON and mtsNFVS on computing an negative feedback vertex set of the SIPC model <https://doi.org/10.7554/eLife.72626>.
    Note that the result of mtsNFVS is not deterministic.
    """
    file_path = os.getcwd() + "/SIPC.bnet"
    bn_real = BooleanNetwork.from_file(file_path)

    nfvs_aeon = feedback_vertex_set(bn_real, parity='negative')
    nfvs_mtsNFVS = find_minimum_NFVS(bn_real)

    assert len(nfvs_mtsNFVS) <= 13 # the result of mtsNFVS is 13
    assert len(nfvs_aeon) <= 25 # the result of AEON is 26
    

