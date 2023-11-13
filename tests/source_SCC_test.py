import sys
sys.path.append("../nfvsmotifs/_sd_algorithms")
sys.path.append("..")


from expand_source_SCCs import *

from nfvsmotifs.SuccessionDiagram import SuccessionDiagram

def test_find_source_nodes():

    bn = BooleanNetwork.from_bnet(
    """targets,factors
    constant1_1, (constant1_1 | !constant1_1)
    constant1_0, (constant1_0 & !constant1_0)
    constant2_1, true
    constant2_0, false
    source, source
    oscillator, !oscillator
    source_after_perc, source_after_perc & constant1_1
    after_perc_0, after_perc_0 & constant1_0""")

    source_nodes = find_source_nodes(bn)

    assert source_nodes == ["source"]

    perc_space, _ = percolate_space(bn, {}, strict_percolation=False)
    perc_bn = percolate_network(bn, perc_space)

    source_nodes = find_source_nodes(perc_bn)

    assert source_nodes == ["source", "source_after_perc"]

test_find_source_nodes()

def test_perc_and_remove_constants_from_bn():
    
    bn = BooleanNetwork.from_bnet(
    """targets,factors
    constant1_1, (constant1_1 | !constant1_1)
    constant1_0, (constant1_0 & !constant1_0)
    constant2_1, true
    constant2_0, false
    source, source
    oscillator, !oscillator
    source_after_perc, source_after_perc & constant1_1
    after_perc_0, after_perc_0 & constant1_0""")

    clean_bnet, clean_bn = perc_and_remove_constants_from_bn(bn, {})

    bn2 = BooleanNetwork.from_bnet(
    """targets,factors
    source, source
    oscillator, !oscillator
    source_after_perc, source_after_perc""")

    clean_bnet2, clean_bn2 = perc_and_remove_constants_from_bn(bn2, {})

    assert clean_bnet == clean_bnet2


test_perc_and_remove_constants_from_bn()

def basic_analysis(bn:BooleanNetwork):
    sd = SuccessionDiagram(bn)
    fully_expanded = expand_source_SCCs(sd)
    assert fully_expanded

    attractor_count = 0
    motif_avoidant_count = 0
    for node in sd.node_ids():
        attr = sd.node_attractor_seeds(node, compute=True)
        attractor_count += len(attr)
        if not sd.node_is_minimal(node):
            motif_avoidant_count += len(attr)

    return bn.num_vars(), len(sd), sd.depth(), attractor_count, motif_avoidant_count, len(sd.minimal_trap_spaces())

def basic_analysis_bfs(bn:BooleanNetwork):
    sd = SuccessionDiagram(bn)
    fully_expanded = sd.expand_bfs()
    assert fully_expanded

    attractor_count = 0
    motif_avoidant_count = 0
    for node in sd.node_ids():
        attr = sd.node_attractor_seeds(node, compute=True)
        attractor_count += len(attr)
        if not sd.node_is_minimal(node):
            motif_avoidant_count += len(attr)

    return bn.num_vars(), len(sd), sd.depth(), attractor_count, motif_avoidant_count, len(sd.minimal_trap_spaces())

def test_expand_source_SCCs():
    bn = BooleanNetwork.from_bnet(
    """targets,factors
    constant1_1, (constant1_1 | !constant1_1)
    constant1_0, (constant1_0 & !constant1_0)
    constant2_1, true
    constant2_0, false
    source, source
    oscillator, !oscillator
    source_after_perc, source_after_perc & constant1_1
    after_perc_0, after_perc_0 & constant1_0""")

    N, size, depth, att, maa, min = basic_analysis(bn)

    assert N == 8
    assert size == 5
    assert depth == 1
    assert att == 4
    assert maa == 0
    assert min == 4

    bn = BooleanNetwork.from_bnet(
    """targets,factors
    source1, source1
    source2, source2
    oscillator, !oscillator""")

    N, size, depth, att, maa, min = basic_analysis(bn)

    assert N == 3
    assert size == 5
    assert depth == 1
    assert att == 4
    assert maa == 0
    assert min == 4

    bn = BooleanNetwork.from_bnet(
    """targets,factors
    constant, true""")

    N, size, depth, att, maa, min = basic_analysis(bn)

    assert N == 1
    assert size == 1
    assert depth == 0
    assert att == 1
    assert maa == 0
    assert min == 1

    bn = BooleanNetwork.from_bnet(
    """targets,factors
    oscillator, !oscillator""")

    N, size, depth, att, maa, min = basic_analysis(bn)

    assert N == 1
    assert size == 1
    assert depth == 0
    assert att == 1
    assert maa == 0
    assert min == 1

    bn = BooleanNetwork.from_bnet(
    """targets,factors
    source1, source1
    source2, source2
    A, B & source1
    B, A
    C, D & source2
    D, C""")

    N, size, depth, att, maa, min = basic_analysis(bn)

    assert N == 6
    assert size == 15
    assert depth == 3
    assert att == 9
    assert maa == 0
    assert min == 9

    bn = BooleanNetwork.from_bnet(
    """targets,factors
    source1, source1
    source2, source2
    A, B & source1
    B, A
    C, D & source2
    D, C""")

    N, size, depth, att, maa, min = basic_analysis(bn)

    assert N == 6
    assert size == 15
    assert depth == 3
    assert att == 9
    assert maa == 0
    assert min == 9

    # motif avoidant attractor
    bn = BooleanNetwork.from_bnet(
    """targets,factors
    A, !A & !B | C
    B, !A & !B | C
    C, A & B""")

    N, size, depth, att, maa, min = basic_analysis(bn)

    assert N == 3
    assert size == 2
    assert depth == 1
    assert att == 2
    assert maa == 1
    assert min == 1

    # motif avoidant attractor2
    bn = BooleanNetwork.from_bnet(
    """targets,factors
    A, !C | (A & B & C)
    B, !A | (A & B & C)
    C, !B | (A & B & C)""")

    N, size, depth, att, maa, min = basic_analysis(bn)

    assert N == 3
    assert size == 2
    assert depth == 1
    assert att == 2
    assert maa == 1
    assert min == 1

    # real bnet
    print("working on 014.bnet")
    bn = BooleanNetwork.from_file("../models/bbm-bnet-inputs-true/014.bnet")

    N, size, depth, att, maa, min = basic_analysis(bn)

    assert att == 2
    assert maa == 0
    assert min == 2

    # real bnet
    print("working on 177.bnet")
    bn = BooleanNetwork.from_file("../models/bbm-bnet-inputs-true/177.bnet")

    N, size, depth, att, maa, min = basic_analysis(bn)

    assert att == 2
    assert maa == 0
    assert min == 2



test_expand_source_SCCs()

def test_find_scc_sd():
    bnet = """targets,factors
A, B
B, A | A &  C"""

    scc_sd = find_scc_sd(bnet, ["A", "B"])

    assert scc_sd.G.nodes[0]["space"] == {}
    assert scc_sd.G.nodes[1]["space"] == {"A":0, "B":0}
    assert scc_sd.G.nodes[2]["space"] == {"A":1, "B":1}


test_find_scc_sd()