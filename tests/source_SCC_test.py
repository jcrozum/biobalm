from biodivine_aeon import BooleanNetwork

from nfvsmotifs._sd_algorithms.expand_source_SCCs import (
    expand_source_SCCs,
    find_scc_sd,
    find_source_nodes,
    perc_and_remove_constants_from_bn,
)
from nfvsmotifs.space_utils import percolate_network, percolate_space
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
    after_perc_0, after_perc_0 & constant1_0"""
    )

    source_nodes = find_source_nodes(bn)

    assert source_nodes == ["source"]

    perc_space, _ = percolate_space(bn, {}, strict_percolation=False)
    perc_bn = percolate_network(bn, perc_space)

    source_nodes = find_source_nodes(perc_bn)

    assert source_nodes == ["source", "source_after_perc"]


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
    after_perc_0, after_perc_0 & constant1_0"""
    )

    clean_bnet, _ = perc_and_remove_constants_from_bn(bn, {})

    bn2 = BooleanNetwork.from_bnet(
        """targets,factors
    source, source
    oscillator, !oscillator
    source_after_perc, source_after_perc"""
    )

    clean_bnet2, _ = perc_and_remove_constants_from_bn(bn2, {})

    assert clean_bnet == clean_bnet2


def test_find_scc_sd():
    bnet = """targets,factors
A, B
B, A | A &  C"""

    scc_sd, _ = find_scc_sd(
        bnet, ["A", "B"], expander=SuccessionDiagram.expand_bfs, check_maa=True
    )

    assert scc_sd.G.nodes[0]["space"] == {}
    assert scc_sd.G.nodes[1]["space"] == {"A": 0, "B": 0}
    assert scc_sd.G.nodes[2]["space"] == {"A": 1, "B": 1}


def expansion(bn: BooleanNetwork):
    sd = SuccessionDiagram(bn)
    fully_expanded = expand_source_SCCs(sd, check_maa=False)
    assert fully_expanded

    return bn.num_vars(), len(sd), sd.depth(), len(sd.minimal_trap_spaces())


def test_expansion():
    # motif avoidant attractor
    bn = BooleanNetwork.from_bnet(
        """targets,factors
    A, !A & !B | C
    B, !A & !B | C
    C, A & B"""
    )
    n, size, depth, min = expansion(bn)
    assert n == 3
    assert size == 2
    assert depth == 1
    assert min == 1

    # motif avoidant attractor2
    bn = BooleanNetwork.from_bnet(
        """targets,factors
    A, !C | (A & B & C)
    B, !A | (A & B & C)
    C, !B | (A & B & C)"""
    )
    n, size, depth, min = expansion(bn)
    assert n == 3
    assert size == 2
    assert depth == 1
    assert min == 1

    # motif avoidant attractors
    bn = BooleanNetwork.from_bnet(
        """targets,factors
    A, !A & !B | C
    B, !A & !B | C
    C, A & B
    X, !Z | (X & Y & Z)
    Y, !X | (X & Y & Z)
    Z, !Y | (X & Y & Z)"""
    )
    n, size, depth, min = expansion(bn)
    assert n == 6
    assert size == 3
    assert depth == 2
    assert min == 1

    # real bnet
    bn = BooleanNetwork.from_file("models/bbm-bnet-inputs-true/002.bnet")
    _, _, _, min = expansion(bn)
    assert min == 24


def attractor_search(bn: BooleanNetwork):
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

    return (
        bn.num_vars(),
        len(sd),
        sd.depth(),
        attractor_count,
        motif_avoidant_count,
        len(sd.minimal_trap_spaces()),
    )


def test_attractor_search():
    bn = BooleanNetwork.from_bnet(
        """targets,factors
    constant1_1, (constant1_1 | !constant1_1)
    constant1_0, (constant1_0 & !constant1_0)
    constant2_1, true
    constant2_0, false
    source, source
    oscillator, !oscillator
    source_after_perc, source_after_perc & constant1_1
    after_perc_0, after_perc_0 & constant1_0"""
    )
    N, size, depth, att, maa, min = attractor_search(bn)
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
    oscillator, !oscillator"""
    )
    n, size, depth, att, maa, min = attractor_search(bn)
    assert n == 3
    assert size == 5
    assert depth == 1
    assert att == 4
    assert maa == 0
    assert min == 4

    bn = BooleanNetwork.from_bnet(
        """targets,factors
    constant, true"""
    )
    n, size, depth, att, maa, min = attractor_search(bn)
    assert n == 1
    assert size == 1
    assert depth == 0
    assert att == 1
    assert maa == 0
    assert min == 1

    bn = BooleanNetwork.from_bnet(
        """targets,factors
    oscillator, !oscillator"""
    )
    n, size, depth, att, maa, min = attractor_search(bn)
    assert n == 1
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
    D, C"""
    )
    n, size, depth, att, maa, min = attractor_search(bn)
    assert n == 6
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
    D, C"""
    )
    n, size, depth, att, maa, min = attractor_search(bn)
    assert n == 6
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
    C, A & B"""
    )
    n, size, depth, att, maa, min = attractor_search(bn)
    assert n == 3
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
    C, !B | (A & B & C)"""
    )
    n, size, depth, att, maa, min = attractor_search(bn)
    assert n == 3
    assert size == 2
    assert depth == 1
    assert att == 2
    assert maa == 1
    assert min == 1

    # motif avoidant attractors
    bn = BooleanNetwork.from_bnet(
        """targets,factors
    A, !A & !B | C
    B, !A & !B | C
    C, A & B
    X, !Z | (X & Y & Z)
    Y, !X | (X & Y & Z)
    Z, !Y | (X & Y & Z)"""
    )
    n, size, depth, att, maa, min = attractor_search(bn)
    assert n == 6
    assert size == 4
    assert depth == 2
    assert att == 4
    assert maa == 3
    assert min == 1

    # complicated combination
    bn = BooleanNetwork.from_bnet(
        """targets,factors
    source, source
    X, Y & source
    Y, X
    A, !A & !B | C
    B, !A & !B | C
    C, A & B & source"""
    )
    n, size, depth, att, maa, min = attractor_search(bn)
    assert n == 6
    assert size == 7
    assert depth == 3
    assert att == 5
    assert maa == 2
    assert min == 3

    # 3 levels
    bn = BooleanNetwork.from_bnet(
        """targets,factors
    X1, Y1
    Y1, X1
    V1, W1
    W1, V1
    A, !A & !B | C
    B, !A & !B | C
    C, A & B
    X2, Y2 & X1
    Y2, X2
    V2, W2 & X1
    W2, V2
    X3, Y3 & X2
    Y3, X3"""
    )
    n, size, depth, att, maa, min = attractor_search(bn)
    assert n == 13
    assert size == 41
    assert depth == 6
    assert att == 28
    assert maa == 14
    assert min == 14

    # real bnet
    bn = BooleanNetwork.from_file("models/bbm-bnet-inputs-true/014.bnet")
    n, size, depth, att, maa, min = attractor_search(bn)
    assert att == 2
    assert maa == 0
    assert min == 2

    # real bnet
    bn = BooleanNetwork.from_file("models/bbm-bnet-inputs-true/177.bnet")
    n, size, depth, att, maa, min = attractor_search(bn)
    assert att == 2
    assert maa == 0
    assert min == 2
