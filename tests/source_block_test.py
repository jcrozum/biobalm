from biodivine_aeon import BooleanNetwork

from biobalm import SuccessionDiagram


def expansion(bn: BooleanNetwork):
    sd = SuccessionDiagram(bn)
    fully_expanded = sd.expand_block(find_motif_avoidant_attractors=False)
    assert fully_expanded

    return (
        bn.variable_count(),
        len(sd),
        len(list(sd.expanded_ids())),
        sd.depth(),
        len(sd.minimal_trap_spaces()),
    )


def test_expansion():
    # motif avoidant attractor
    bn = BooleanNetwork.from_bnet(
        """targets,factors
    A, !A & !B | C
    B, !A & !B | C
    C, A & B"""
    )
    n, size, e_size, depth, min = expansion(bn)
    assert n == 3
    assert size == 2
    assert e_size == 2
    assert depth == 1
    assert min == 1

    # motif avoidant attractor2
    bn = BooleanNetwork.from_bnet(
        """targets,factors
    A, !C | (A & B & C)
    B, !A | (A & B & C)
    C, !B | (A & B & C)"""
    )
    n, size, e_size, depth, min = expansion(bn)
    assert n == 3
    assert size == 2
    assert e_size == 2
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
    n, size, e_size, depth, min = expansion(bn)
    assert n == 6
    assert size == 4
    assert e_size == 3
    assert depth == 2
    assert min == 1

    # real bnet
    bn = BooleanNetwork.from_file("models/bbm-bnet-inputs-true/002.bnet")
    _, _, _, _, min = expansion(bn)
    assert min == 24


def attractor_search(bn: BooleanNetwork):
    sd = SuccessionDiagram(bn)
    fully_expanded = sd.expand_block(find_motif_avoidant_attractors=True)
    assert fully_expanded

    attractor_count = 0
    motif_avoidant_count = 0
    for node in sd.expanded_ids():
        attr = sd.node_attractor_seeds(node, compute=True)
        attractor_count += len(attr)
        if not sd.node_is_minimal(node):
            motif_avoidant_count += len(attr)

    return (
        bn.variable_count(),
        len(sd),
        len(list(sd.expanded_ids())),
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
    N, size, e_size, depth, att, maa, min = attractor_search(bn)
    assert N == 8
    assert size == 5
    assert e_size == 5
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
    n, size, e_size, depth, att, maa, min = attractor_search(bn)
    assert n == 3
    assert size == 5
    assert e_size == 5
    assert depth == 1
    assert att == 4
    assert maa == 0
    assert min == 4

    bn = BooleanNetwork.from_bnet(
        """targets,factors
    constant, true"""
    )
    n, size, e_size, depth, att, maa, min = attractor_search(bn)
    assert n == 1
    assert size == 1
    assert e_size == 1
    assert depth == 0
    assert att == 1
    assert maa == 0
    assert min == 1

    bn = BooleanNetwork.from_bnet(
        """targets,factors
    oscillator, !oscillator"""
    )
    n, size, e_size, depth, att, maa, min = attractor_search(bn)
    assert n == 1
    assert size == 1
    assert e_size == 1
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
    n, size, e_size, depth, att, maa, min = attractor_search(bn)
    assert n == 6
    assert size == 17
    assert e_size == 15
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
    n, size, e_size, depth, att, maa, min = attractor_search(bn)
    assert n == 6
    assert size == 17
    assert e_size == 15
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
    n, size, e_size, depth, att, maa, min = attractor_search(bn)
    assert n == 3
    assert size == 2
    assert e_size == 2
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
    n, size, e_size, depth, att, maa, min = attractor_search(bn)
    assert n == 3
    assert size == 2
    assert e_size == 2
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
    n, size, e_size, depth, att, maa, min = attractor_search(bn)
    assert n == 6
    assert size == 4
    assert e_size == 4
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
    n, size, e_size, depth, att, maa, min = attractor_search(bn)
    assert n == 6
    assert size == 8
    assert e_size == 7
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
    n, size, e_size, depth, att, maa, min = attractor_search(bn)
    assert n == 13
    assert size == 75
    assert e_size == 41
    assert depth == 6
    assert att == 28
    assert maa == 14
    assert min == 14

    # real bnet
    bn = BooleanNetwork.from_file("models/bbm-bnet-inputs-true/014.bnet")
    n, size, e_size, depth, att, maa, min = attractor_search(bn)
    assert att == 2
    assert maa == 0
    assert min == 2

    # real bnet
    bn = BooleanNetwork.from_file("models/bbm-bnet-inputs-true/177.bnet")
    n, size, e_size, depth, att, maa, min = attractor_search(bn)
    assert att == 2
    assert maa == 0
    assert min == 2

    # interesting example uncovered during testing (random-nk2/n20_29.bnet)

    bn = BooleanNetwork.from_bnet(
        """targets, factors
    n0, (n12 & !n5) | (n12 & n5)
    n1, (!n1 & !n15) | (n1 & !n15) | (n1 & n15)
    n2, (n3 & !n15) | (n3 & n15)
    n3, (!n17 & n3) | (n17 & !n3)
    n4, (n0 & !n18) | (n0 & n18)
    n5, (!n13 & !n16) | (!n13 & n16) | (n13 & n16)
    n6, true
    n7, (!n10 & !n1) | (n10 & !n1) | (n10 & n1)
    n8, true
    n9, (!n0 & !n19) | (!n0 & n19) | (n0 & !n19)
    n10, (!n2 & !n17) | (!n2 & n17)
    n11, (n16 & !n1)
    n12, (!n3 & n12) | (n3 & n12)
    n13, (!n8 & n3) | (n8 & !n3)
    n14, (!n3 & !n15) | (n3 & !n15)
    n15, (!n14 & !n15) | (!n14 & n15)
    n16, (!n15 & !n11) | (n15 & !n11)
    n17, (n2 & !n11) | (n2 & n11)
    n18, (!n7 & !n18) | (!n7 & n18)
    n19, (!n5 & !n12) | (!n5 & n12)"""
    )
    n, size, e_size, depth, att, maa, min = attractor_search(bn)
    assert size == 21
    assert e_size == 13
    assert min == 6
    assert att == 6
    assert maa == 0


def test_isomorph():
    path = "models/bbm-bnet-inputs-true/005.bnet"
    bn = BooleanNetwork.from_file(path)

    sd_bfs = SuccessionDiagram(bn)
    sd_bfs.expand_bfs()

    sd_scc = SuccessionDiagram(bn)
    sd_scc.expand_scc()

    assert [sd_bfs.node_data(id)["space"] for id in sd_bfs.node_ids()] == [
        sd_scc.node_data(id)["space"] for id in sd_scc.node_ids()
    ]

    assert sd_scc.is_isomorphic(sd_bfs)

    edge_motifs_bfs = set(
        str(sorted(sd_bfs.edge_stable_motif(x, y).items()))  # type: ignore
        for (x, y) in sd_bfs.dag.edges  # type: ignore
    )  # type: ignore
    edge_motifs_scc = set(
        str(sorted(sd_scc.edge_stable_motif(x, y).items()))  # type: ignore
        for (x, y) in sd_scc.dag.edges  # type: ignore
    )  # type: ignore

    assert edge_motifs_bfs == edge_motifs_scc
