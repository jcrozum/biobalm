from biodivine_aeon import BooleanNetwork  # type:ignore
from networkx import DiGraph  # type:ignore

from nfvsmotifs.interaction_graph_utils import (
    feedback_vertex_set,
    find_minimum_NFVS,
    independent_cycles,
    infer_signed_interaction_graph,
)


def test_ig_inference():
    bn = BooleanNetwork.from_bnet(
        """
        # Just a normal function.
        b, a | !b
        # Contradiciton on `a` - the regulation should not appear in the result
        # Also, non-monotonic dependence on b and c.
        a, (a & !a) | (b <=> c)
        c, c
    """
    )
    ig = infer_signed_interaction_graph(bn)

    edges = {edge: ig.get_edge_data(edge[0], edge[1])["sign"] for edge in ig.edges}  # type: ignore
    assert len(edges) == 5  # type: ignore
    assert edges[("a", "b")] == "+"
    assert edges[("b", "b")] == "-"
    assert edges[("b", "a")] == "?"
    assert edges[("c", "a")] == "?"
    assert edges[("c", "c")] == "+"
    assert ("a", "a") not in edges


# There should be a negative cycle between b_1 and b_2,
# a positive cycle between d_1 and d_2, and a negative cycle
# between d_1, d_2, and d_3. Other nodes are not on cycles
# except for e, which has a positive self-loop.
CYCLES_BN = BooleanNetwork.from_aeon(
    """
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
    """
)

CYCLES_DIGRAPH = DiGraph()
CYCLES_DIGRAPH.add_nodes_from(["a", "b_1", "b_2", "c", "d_1", "d_2", "d_3", "e"])  # type: ignore
CYCLES_DIGRAPH.add_edges_from(  # type: ignore
    [
        ("a", "c", {"sign": "+"}),
        ("b_1", "b_2", {"sign": "+"}),
        ("b_2", "b_1", {"sign": "-"}),
        ("b_2", "c", {"sign": "+"}),
        ("c", "d_2", {"sign": "+"}),
        ("c", "e", {"sign": "+"}),
        ("d_1", "d_3", {"sign": "+"}),
        ("d_3", "d_2", {"sign": "-"}),
        ("d_2", "d_1", {"sign": "+"}),
        ("d_1", "d_2", {"sign": "+"}),
        ("e", "e", {"sign": "+"}),
    ]
)


def test_fvs():
    fvs = feedback_vertex_set(CYCLES_BN)
    nfvs = feedback_vertex_set(CYCLES_BN, parity="negative")
    pfvs = feedback_vertex_set(CYCLES_BN, parity="positive")

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
    dg_nfvs = feedback_vertex_set(CYCLES_DIGRAPH, parity="negative")
    dg_pfvs = feedback_vertex_set(CYCLES_DIGRAPH, parity="positive")

    assert fvs == dg_fvs
    assert nfvs == dg_nfvs
    assert pfvs == dg_pfvs


def test_subgraph_fvs():
    # We only keep the two cycles consisting of "d_*". The "b_*" cycle
    # and "e" self-loop are not considered.
    fvs = feedback_vertex_set(CYCLES_BN, subgraph=["a", "b_1", "d_1", "d_2", "d_3"])
    pfvs = feedback_vertex_set(
        CYCLES_BN, parity="positive", subgraph=["a", "b_1", "d_1", "d_2", "d_3"]
    )
    nfvs = feedback_vertex_set(
        CYCLES_BN, parity="negative", subgraph=["a", "b_1", "d_1", "d_2", "d_3"]
    )

    assert len(fvs) == 1
    assert len(pfvs) == 1
    assert len(nfvs) == 1
    assert ("d_1" in fvs) or ("d_2" in fvs)
    assert ("d_1" in nfvs) or ("d_2" in nfvs) or ("d_3" in nfvs)
    assert ("d_1" in pfvs) or ("d_2" in pfvs)


def test_ic():
    ic = independent_cycles(CYCLES_BN)
    n_ic = independent_cycles(CYCLES_BN, parity="negative")
    p_ic = independent_cycles(CYCLES_BN, parity="positive")

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
    dg_ic = independent_cycles(CYCLES_DIGRAPH)  # type: ignore
    dg_n_ic = independent_cycles(CYCLES_DIGRAPH, parity="negative")  # type: ignore
    dg_p_ic = independent_cycles(CYCLES_DIGRAPH, parity="positive")  # type: ignore

    print(ic)
    print(dg_ic)

    assert [set(x) for x in ic] == [set(x) for x in dg_ic]
    assert [set(x) for x in n_ic] == [set(x) for x in dg_n_ic]
    assert [set(x) for x in p_ic] == [set(x) for x in dg_p_ic]


def test_subgraph_ic():
    # We only keep the two cycles consisting of "d_*". The "b_*" cycle
    # and "e" self-loop are not considered.
    ic = independent_cycles(CYCLES_BN, subgraph=["a", "b_1", "d_1", "d_2", "d_3"])
    p_ic = independent_cycles(
        CYCLES_BN, parity="positive", subgraph=["a", "b_1", "d_1", "d_2", "d_3"]
    )
    n_ic = independent_cycles(
        CYCLES_BN, parity="negative", subgraph=["a", "b_1", "d_1", "d_2", "d_3"]
    )

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
    bn_real = BooleanNetwork.from_bnet(
        """
        ABL1, (ATM & !RB1)
        ACVR1, BMPR2
        ADAM17, ERK_f
        AKT1S1, !AKT_f
        AKT_f, (!ATM & ((!PDPK1 & (!PHLPP1 & (!PPP1CA & mTORC2_c))) | (PDPK1 & (!PHLPP1 & !PPP1CA))))
        AP1_c, ((!FOS & SMAD4) | FOS)
        APAF1, CYCS
        APC, ((!AXIN1 & (!DVL_f & (GSK3_f & !PRKACA))) | (AXIN1 & (!DVL_f & !PRKACA)))
        ARHGAP24, ROCK1
        ATF2, ((!ERK_f & MAPK14) | ERK_f)
        ATM, ((!ATR & BRCA1) | ATR)
        ATR, ((!ABL1 & CDKN2A) | ABL1)
        AURKA, PAK1
        AURKB, (!ATM & CHEK1)
        AXIN1, (!AURKA & (GSK3_f & (!LRP_f & (!PPM1A & !PPP1CA))))
        BAD, (!AKT_f & (!PAK1 & (!RAF_f & !RSK_f)))
        BAK1, (!BCL2 & ((!BID & (!MCL1 & TP53)) | (BID & !MCL1)))
        BAX, ((!BID & (!MCL1 & TP53)) | (BID & !MCL1))
        BBC3, TP53
        BCL2, (!BAD & (!BBC3 & (!BID & (CREB1 & !TP53))))
        BID, ((!CASP8 & CSNK1A1) | CASP8)
        BIRC_f, (!AURKB & (!DIABLO & (STAT3 & !TP53)))
        BMPR2, (!SMURF1 & !SMURF2)
        BRCA1, ((!AKT_f & (!CCND1 & CHEK2)) | (AKT_f & !CCND1))
        BTRC, ((!AXIN1 & (GSK3_f & !LRP_f)) | (AXIN1 & !LRP_f))
        BUB1, ATM
        CASP3, (!BIRC_f & ((!CASP8 & (CASP9 & !XIAP)) | (CASP8 & !XIAP)))
        CASP8, (!CFLAR & (DD_f & !SRC))
        CASP9, (!BIRC_f & (PPP1CA & !XIAP))
        CBPp300_c, (EP300 & !TP53)
        CCNB1, !CDKN1A
        CCND1, (!CDKN2A & (!CHUK & (!GSK3_f & ((!RSK_f & STAT3) | RSK_f))))
        CCNE1, ((!CDC25A & (!CDKN1A & (!CDKN1B & (E2F1 & !TGFB1)))) | (CDC25A & (!CDKN1A & (!CDKN1B & !TGFB1))))
        CDC25A, (!CHEK1 & (!CHEK2 & (!CSNK1A1 & (!GSK3_f & (MYC & !SMAD3)))))
        CDC42, (!ARHGAP24 & SRC)
        CDH1, (!SNAI_f & !TWIST1)
        CDKN1A, (!AKT_f & ((!BRCA1 & (!MDM2 & (!MYC & (!SKP2 & TP53)))) | (BRCA1 & (!MDM2 & (!MYC & !SKP2)))))
        CDKN1B, (!AKT_f & (!MYC & TGFB1))
        CDKN2A, !MYC
        CFLAR, ((!AKT_f & (!ITCH & MAPK14)) | (AKT_f & !ITCH))
        CFL_f, ((!LIMK1 & LIMK2) | LIMK1)
        CHEK1, (!AKT_f & ATR)
        CHEK2, ((!ATR & PLK1) | ATR)
        CHUK, ((!AKT_f & TRAF6) | AKT_f)
        CREB1, ((!AKT_f & (!ATM & (!ATR & SMAD2))) | (AKT_f & (!ATM & !ATR)))
        CREBBP, CHUK
        CSK, PRKACA
        CSNK1A1, !LRP_f
        CSNK1_f, ((!AXIN1 & DVL_f) | AXIN1)
        CTNNB1, (!APC & (!BTRC & ((!CHUK & (!CSNK1A1 & (!CSNK1_f & YAP_TAZ))) | (CHUK & (!CSNK1A1 & !CSNK1_f)))))
        CYCS, ((!BAK1 & (BAX & !BCL2)) | (BAK1 & !BCL2))
        DAAM1, DVL_f
        DD_f, (AURKA & (!CSNK1A1 & !MAP2K7))
        DIABLO, ((!BAK1 & (BAX & (!BCL2 & !BIRC_f))) | (BAK1 & (!BCL2 & !BIRC_f)))
        DKK_f, (!MYC & TCF7_f)
        DUSP1, ((!ERK_f & (MSK_f & !SKP2)) | (ERK_f & !SKP2))
        DUSP6, ((!ERK_f & mTORC1_c) | ERK_f)
        DVL_f, ((!FZD_f & (!ITCH & (SMAD1 & !YAP_TAZ))) | (FZD_f & (!ITCH & !YAP_TAZ)))
        E2F1, (!HES1 & !RB1)
        EGR1, (NFKB_f & !TCF7_f)
        EP300, (AKT_f & (!PRKCD & !SKI))
        ERK_f, (!DUSP6 & (MEK_f & !PPP1CA))
        FOS, ((!ERK_f & SRF) | ERK_f)
        FOXO_f, (!AKT_f & (!CREB1 & (!CSNK1A1 & (DUSP6 & (!IKBKB & !NLK)))))
        FZD_f, !SFRP1
        GAB_f, (!ERK_f & GRB2)
        GLI_f, (!CSNK1A1 & (!GSK3_f & (!PRKACA & !SUFU)))
        GRB2, ((!RTPK_f & SHC1) | RTPK_f)
        GSK3_f, (!AKT_f & (CSNK1A1 & (!DVL_f & (!ERK_f & (!MAPK14 & (!RSK_f & !S6K_f))))))
        HES1, ((!NOTCH1 & REL_f) | NOTCH1)
        IKBKB, (DD_f & (!PLK1 & (!PPM1A & !TP53)))
        ILK, (PAK1 & !TWIST1)
        ILR_f, ((!AP1_c & LIF) | AP1_c)
        IRAK1, (ILR_f & !SOCS1)
        IRS1, (!ERK_f & (!IKBKB & !S6K_f))
        ITCH, JNK_f
        JAK_f, (ILR_f & (!PTPN6 & !SOCS1))
        JNK_f, (!DUSP1 & ((!MAP2K4 & MAP2K7) | MAP2K4))
        JUN, (!GSK3_f & JNK_f)
        KRAS, ((!PTPN11 & SOS1) | PTPN11)
        LATS_f, ((!AURKA & MOB1_f) | AURKA)
        LEF1, (!CSNK1_f & CTNNB1)
        LIF, RAF_f
        LIMK1, ((!RAC_f & ROCK1) | RAC_f)
        LIMK2, (!PRKCD & ROCK1)
        LRP_f, (!DKK_f & ((!ERK_f & FZD_f) | ERK_f))
        MAP2K3, ((!MAP3K5 & MAP3K7) | MAP3K5)
        MAP2K4, ((!MAP3K11 & MAP3K4) | MAP3K11)
        MAP2K7, ((!MAP3K7 & MAPK8IP3) | MAP3K7)
        MAP3K11, RAC_f
        MAP3K4, RAC_f
        MAP3K5, !AKT_f
        MAP3K7, ((!TAB_f & TRAF6) | TAB_f)
        MAP3K8, IKBKB
        MAPK14, (!DUSP1 & ((!MAP2K3 & MAP2K4) | MAP2K3))
        MAPK8IP3, ROCK1
        MAPKAPK2, MAPK14
        MCL1, (!BBC3 & ((!ERK_f & (!GSK3_f & JNK_f)) | (ERK_f & !GSK3_f)))
        MDM2, (!ABL1 & ((!AKT_f & (!ATM & (!ATR & (!CDKN2A & (!CSNK1_f & (!S6K_f & TP53)))))) | (AKT_f & (!ATM & (!ATR & (!CDKN2A & (!CSNK1_f & !S6K_f)))))))
        MEK_f, (!ERK_f & RAF_f)
        MMP_f, ((!LEF1 & STAT3) | LEF1)
        MOB1_f, STK_f
        MSK_f, ((!ERK_f & MAPK14) | ERK_f)
        MYC, (!GSK3_f & ((!PLK1 & STAT3) | PLK1))
        NFKB_f, (!CHEK1 & ((!IKBKB & MSK_f) | IKBKB))
        NLK, MAP3K7
        NOTCH1, ADAM17
        PAK1, ((!CDC42 & RAC_f) | CDC42)
        PARD6A, ((!TGFBR1 & TGFBR2) | TGFBR1)
        PDPK1, (PIK3CA & !PTEN)
        PHLPP1, !GSK3_f
        PIAS1, MAPKAPK2
        PIK3CA, ((!GAB_f & KRAS) | GAB_f)
        PLCG1, SYK
        PLK1, ((!MAPKAPK2 & PDPK1) | MAPKAPK2)
        PPM1A, PTEN
        PPP1CA, (!RTPK_f & SMAD7)
        PRKACA, ((!FOS & NFKB_f) | FOS)
        PRKCA, ((!PARD6A & (!PHLPP1 & PLCG1)) | (PARD6A & !PHLPP1))
        PRKCD, ((!CASP3 & PDPK1) | CASP3)
        PRKDC, ((!ATM & ATR) | ATM)
        PTCH1, GLI_f
        PTEN, (!CBPp300_c & (!GSK3_f & (ROCK1 & !SRC)))
        PTPN11, GAB_f
        PTPN6, (!PRKCA & SRC)
        RAC_f, (!ARHGAP24 & ((!PIK3CA & TIAM1) | PIK3CA))
        RAF_f, (!AKT_f & (!ERK_f & (KRAS & !RHEB)))
        RB1, (!CCND1 & (!CCNE1 & ((!CHEK1 & CHEK2) | CHEK1)))
        RBPJ, !HES1
        REL_f, ((!IKBKB & (PRKCA & !STAT1)) | (IKBKB & !STAT1))
        RHEB, !TSC_f
        RHOA, (DAAM1 & (!PARD6A & (!RAC_f & (!RND3 & !SMURF1))))
        RND3, ROCK1
        ROCK1, ((!CASP3 & RHOA) | CASP3)
        RSK_f, ((!ERK_f & PDPK1) | ERK_f)
        RTPK_f, ((!FOXO_f & (!MAPK14 & (!MEK_f & MMP_f))) | (FOXO_f & (!MAPK14 & !MEK_f)))
        S6K_f, ((!PDPK1 & (!PHLPP1 & mTORC1_c)) | (PDPK1 & !PHLPP1))
        SFRP1, !MYC
        SHC1, (!PTEN & ((!SRC & TGFBR1) | SRC))
        SKI, !AKT_f
        SKP2, ((!EP300 & ERK_f) | EP300)
        SMAD1, ((!ACVR1 & (!ERK_f & (!GSK3_f & (!PPM1A & (!SKI & (!SMAD6 & (!SMURF1 & YAP_TAZ))))))) | (ACVR1 & (!ERK_f & (!GSK3_f & (!PPM1A & (!SKI & (!SMAD6 & !SMURF1)))))))
        SMAD2, ((!ACVR1 & (!PPM1A & (!SKI & (!SMURF2 & TGFBR1)))) | (ACVR1 & (!PPM1A & (!SKI & !SMURF2))))
        SMAD3, ((!ACVR1 & (!AKT_f & (!ERK_f & (!GSK3_f & (!PPM1A & (!SKI & (!SMAD6 & (!SMAD7 & TGFBR1)))))))) | (ACVR1 & (!AKT_f & (!ERK_f & (!GSK3_f & (!PPM1A & (!SKI & (!SMAD6 & !SMAD7))))))))
        SMAD4, (!SKI & ((!SMAD2 & (SMAD5 & (!SMAD6 & (!SMAD7 & !SMURF1)))) | (SMAD2 & (!SMAD6 & (!SMAD7 & !SMURF1)))))
        SMAD5, (ACVR1 & (!SKI & !SMURF2))
        SMAD6, ((!SMAD2 & SMAD4) | SMAD2)
        SMAD7, (!AXIN1 & (!ITCH & ((!SMAD2 & (SMAD4 & !SMURF2)) | (SMAD2 & !SMURF2))))
        SMO, (CSNK1A1 & !PTCH1)
        SMURF1, SMAD7
        SMURF2, SMAD7
        SNAI_f, (!GSK3_f & ((!LATS_f & PAK1) | LATS_f))
        SOCS1, STAT1
        SOS1, (!ERK_f & PLCG1)
        SRC, (!CSK & ((!PRKACA & RTPK_f) | PRKACA))
        SRF, ((!CFL_f & RSK_f) | CFL_f)
        STAT1, ((!IKBKB & (!PIAS1 & SRC)) | (IKBKB & !PIAS1))
        STAT3, ((!IRAK1 & (!PPP1CA & SRC)) | (IRAK1 & !PPP1CA))
        STK_f, (!AKT_f & (PHLPP1 & !RAF_f))
        SUFU, !SMO
        SYK, (!CHEK1 & ILR_f)
        TAB_f, (!MAPK14 & TRAF6)
        TCF7_f, (CTNNB1 & !NLK)
        TGFB1, ((!FOS & NFKB_f) | FOS)
        TGFBR1, (!SMAD6 & (!SMAD7 & (!SMURF1 & (!SMURF2 & TGFBR2))))
        TGFBR2, (!SMURF1 & (!SMURF2 & TGFB1))
        TIAM1, !ROCK1
        TP53, (!AURKB & (EP300 & !MDM2))
        TRAF6, ((!IRAK1 & TGFBR1) | IRAK1)
        TSC_f, (!AKT_f & (!ERK_f & (GSK3_f & (!IKBKB & !RSK_f))))
        TWIST1, ((!ERK_f & JNK_f) | ERK_f)
        VAV1, SYK
        XIAP, ((!AKT_f & (BIRC_f & !DIABLO)) | (AKT_f & !DIABLO))
        YAP_TAZ, (!BTRC & (!CSNK1_f & !LATS_f))
        mTORC1_c, (!AKT1S1 & ((!RHEB & RSK_f) | RHEB))
        mTORC2_c, ((!PIK3CA & (!S6K_f & TSC_f)) | (PIK3CA & !S6K_f))
    """
    )

    nfvs_mtsNFVS = find_minimum_NFVS(bn_real)

    assert len(nfvs_mtsNFVS) <= 19  # the result of mtsNFVS is 19

    for _i in range(10):
        nfvs = find_minimum_NFVS(bn_real)
        assert nfvs == nfvs_mtsNFVS


def test_fvs_accuracy_SIPC():
    """
    Compare results of AEON and mtsNFVS on computing an negative feedback vertex set of the SIPC model <https://doi.org/10.7554/eLife.72626>.
    Note that the result of mtsNFVS is not deterministic.
    """
    bn_real = BooleanNetwork.from_bnet(
        """
        AKT, ((!HSPs&(PIP3&!PTCH1))|(HSPs&!PTCH1))
        AMPK, ((!AMP_ATP&(!ATM&(!ATR&(!EGFR&(!FGFR3&HIF1)))))|((!AMP_ATP&(!ATM&(ATR&(!EGFR&!FGFR3))))|((!AMP_ATP&(ATM&(!EGFR&!FGFR3)))|(AMP_ATP&(!EGFR&!FGFR3)))))
        AMP_ATP, !Nutrients
        APAF1, ((!AKT&(!Bak&(!BAX&(!BCL2&(!Bcl_XL&(!Caspase8&(!HSPs&p53)))))))|((!AKT&(!Bak&(!BAX&(!BCL2&(!Bcl_XL&(!Caspase8&HSPs))))))|((!AKT&(!Bak&(!BAX&(!BCL2&(!Bcl_XL&Caspase8)))))|((!AKT&(!Bak&(BAX&(!BCL2&!Bcl_XL))))|(!AKT&(Bak&(!BCL2&!Bcl_XL)))))))
        AR, ((Androgen&(!EP300&(!EZH2&(!GLI&(!HSPs&(!MDM2&(!NCOA3&(!NCOR1&(!NCOR2&(!PKC&(!PTEN&(!SMAD&NKX3_1))))))))))))|((Androgen&(!EP300&(!EZH2&(!GLI&(!HSPs&(!MDM2&(!NCOA3&(!NCOR1&(!NCOR2&(!PKC&(!PTEN&SMAD)))))))))))|((Androgen&(!EP300&(!EZH2&(!GLI&(!HSPs&(!MDM2&(!NCOA3&(!NCOR1&(!NCOR2&(PKC&!PTEN))))))))))|((Androgen&(!EP300&(!EZH2&(!GLI&(!HSPs&(!MDM2&(NCOA3&(!NCOR1&(!NCOR2&!PTEN)))))))))|((Androgen&(!EP300&(!EZH2&(!GLI&(HSPs&(!MDM2&(!NCOR1&(!NCOR2&!PTEN))))))))|((Androgen&(!EP300&(!EZH2&(GLI&(!MDM2&(!NCOR1&(!NCOR2&!PTEN)))))))|((Androgen&(!EP300&(EZH2&(!MDM2&(!NCOR1&(!NCOR2&!PTEN))))))|(Androgen&(EP300&(!MDM2&(!NCOR1&(!NCOR2&!PTEN))))))))))))
        AR_ERG, (AR&fused_event)
        ATM, DNA_Damage
        ATR, ((!DNA_Damage&p14ARF)|DNA_Damage)
        AXIN1, GSK3
        Acidosis, Acidosis
        Androgen, Androgen
        BAD, (!AKT&(!HIF1&(!p90RSK&(!RAF&!YWHAZ))))
        BAX, ((!BCL2&(!ETS1&(!HIF1&(!JNK&(p53&!YWHAZ)))))|(!BCL2&(!ETS1&(!HIF1&(JNK&!YWHAZ)))))
        BCL2, ((!BAD&(!BAX&(!DAXX&(!HSPs&(!NF_kB&(!p53&RUNX2))))))|((!BAD&(!BAX&(!DAXX&(!HSPs&(NF_kB&!p53)))))|(!BAD&(!BAX&(!DAXX&(HSPs&!p53))))))
        BIRC5, (HSPs&!p53)
        BRCA1, ((!ATM&(!ATR&(!CyclinD&(E2F1&!Caspase3))))|((!ATM&(ATR&(!CyclinD&!Caspase3)))|(ATM&(!CyclinD&!Caspase3))))
        Bak, Caspase8
        Bcl_XL, (AR&(!BAD&!p53))
        CHK1_2, (ATM&(ATR&BRCA1))
        COX4I2, HIF1
        Carcinogen, Carcinogen
        Caspase3, ((!Caspase8&Caspase9)|Caspase8)
        Caspase8, (!cFLAR&FADD)
        Caspase9, ((!AKT&(!BIRC5&(CytoC&(APAF1&PTCH1))))|(!AKT&(BIRC5&(CytoC&APAF1))))
        CyclinB, (E2F1&!p21)
        CyclinD, ((!FOXO&(!GLI&(!GSK3&(MYC_MAX&(!NF_kB&(!p15&(!p21&TCF)))))))|((!FOXO&(!GLI&(!GSK3&(MYC_MAX&(NF_kB&(!p15&!p21))))))|(!FOXO&(GLI&(!GSK3&(MYC_MAX&(!p15&!p21)))))))
        CytoC, ((!AKT&(!Bak&(!BAX&(!BCL2&(!Bcl_XL&(!Caspase8&p53))))))|((!AKT&(!Bak&(!BAX&(!BCL2&(!Bcl_XL&Caspase8)))))|((!AKT&(!Bak&(BAX&(!BCL2&!Bcl_XL))))|(!AKT&(Bak&(!BCL2&!Bcl_XL))))))
        DAXX, (!ATM&(!ATR&!SPOP))
        DNA_Damage, (Carcinogen&!SPOP)
        Dsh, WNT
        E2F1, (AR&(!MXI1&!RB1))
        EGF, EGF
        EGFR, ((!Androgen&(EGF&(!FRS2&!TGFBR)))|(Androgen&(!FRS2&!TGFBR)))
        EP300, AKT
        ERG, ERK
        ERK, ((!MEK1_2&RAF)|MEK1_2)
        ETS1, ((!ERK&(!p53&RTK))|(ERK&!p53))
        EZH2, ((AKT&(E2F1&(!ERG&AR_ERG)))|(AKT&(E2F1&ERG)))
        E_cadherin, ((!beta_catenin&(!ERG&(!NF_kB&(!Slug&(!Snail&AR_ERG)))))|(beta_catenin&(!ERG&(!NF_kB&(!Slug&!Snail)))))
        FADD, ((!MAP3K1_3&(!PTCH1&(!SHH&TNFalpha)))|(!MAP3K1_3&(PTCH1&!SHH)))
        FGF, FGF
        FGFR3, (FGF&(!FRS2&(!PKC&!TGFBR)))
        FOS, ((!AR&(!ERK&ETS1))|((!AR&ERK)|AR))
        FOXA1, FOXO
        FOXO, (!AKT&JNK)
        FRS2, (!ERK&(EGFR&(FGFR3&(!FRS2&!TGFBR))))
        GADD45, ((!p53&SMAD)|p53)
        GLI, ((!SMO&(!SPOP&WNT))|(SMO&!SPOP))
        GSH, ((!MYC_MAX&(!NF_kB&p21))|((!MYC_MAX&NF_kB)|MYC_MAX))
        GSK3, (!AKT&(!Dsh&(MEK1_2&(!mTORC2&!p90RSK))))
        HIF1, ((!FOXO&(!HSPs&(!Hypoxia&(!mTORC2&(!MYC_MAX&(!p53&(!PHDs&(!VHL&mTORC1))))))))|((!FOXO&(!HSPs&(!Hypoxia&(!mTORC2&(MYC_MAX&(!p53&(!PHDs&!VHL)))))))|((!FOXO&(!HSPs&(!Hypoxia&(mTORC2&(!p53&(!PHDs&!VHL))))))|((!FOXO&(!HSPs&(Hypoxia&(!p53&(!PHDs&!VHL)))))|(!FOXO&(HSPs&(!p53&(!PHDs&!VHL))))))))
        HSPs, ((!AKT&(!BRCA1&(!Carcinogen&(!E2F1&(!FOXA1&PKC)))))|((!AKT&(!BRCA1&(!Carcinogen&(!E2F1&FOXA1))))|(!AKT&(!BRCA1&(Carcinogen&!E2F1)))))
        Hypoxia, Hypoxia
        IKK, ((!AKT&(!ETS1&(mTORC2&(!p53&(!PHDs&(!PKC&TAK1))))))|((!AKT&(!ETS1&(mTORC2&(!p53&(!PHDs&PKC)))))|((!AKT&(ETS1&(mTORC2&(!p53&!PHDs))))|((AKT&(!ETS1&(!p53&(!PHDs&(!PKC&TAK1)))))|((AKT&(!ETS1&(!p53&(!PHDs&PKC))))|(AKT&(ETS1&(!p53&!PHDs))))))))
        JNK, ((!ATM&(!ERK&(!GADD45&(!MAP3K1_3&(!p38&TAK1)))))|((!ATM&(!ERK&(!GADD45&(MAP3K1_3&!p38))))|((!ATM&(!ERK&(GADD45&!p38)))|(ATM&(!ERK&!p38)))))
        JUN, ((!AR&(!ETS1&JNK))|((!AR&ETS1)|AR))
        MAP3K1_3, (!FOS&(!JUN&(!JNK&(!p38&RAS))))
        MDM2, ((!AKT&(!ATM&(!ATR&(DAXX&(!p14ARF&p53)))))|(AKT&(!ATM&(!ATR&(DAXX&!p14ARF)))))
        MED12, !GLI
        MEK1_2, ((!MAP3K1_3&(!PDK1&RAF))|((!MAP3K1_3&PDK1)|MAP3K1_3))
        MXI1, HIF1
        MYC, ((!E2F1&(!ERG&(!ERK&(!FOS&(!JUN&(!GLI&(!HIF1&(!NF_kB&(!p38&(!SHH&(!TCF&(AR_ERG&mTORC1))))))))))))|((!E2F1&(!ERG&(!ERK&(!FOS&(!JUN&(!GLI&(!HIF1&(!NF_kB&(!p38&(!SHH&(TCF&mTORC1)))))))))))|((!E2F1&(!ERG&(!ERK&(!FOS&(!JUN&(!GLI&(!HIF1&(!NF_kB&(!p38&(SHH&mTORC1))))))))))|((!E2F1&(!ERG&(!ERK&(!FOS&(!JUN&(!GLI&(!HIF1&(!NF_kB&(p38&(!SHH&(!TCF&AR_ERG)))))))))))|((!E2F1&(!ERG&(!ERK&(!FOS&(!JUN&(!GLI&(!HIF1&(!NF_kB&(p38&(!SHH&TCF))))))))))|((!E2F1&(!ERG&(!ERK&(!FOS&(!JUN&(!GLI&(!HIF1&(!NF_kB&(p38&SHH)))))))))|((!E2F1&(!ERG&(!ERK&(!FOS&(!JUN&(!GLI&(!HIF1&(NF_kB&(!p38&mTORC1)))))))))|((!E2F1&(!ERG&(!ERK&(!FOS&(!JUN&(!GLI&(!HIF1&(NF_kB&p38))))))))|((!E2F1&(!ERG&(!ERK&(!FOS&(!JUN&(GLI&(!HIF1&(!p38&mTORC1))))))))|((!E2F1&(!ERG&(!ERK&(!FOS&(!JUN&(GLI&(!HIF1&p38)))))))|((!E2F1&(!ERG&(!ERK&(!FOS&(JUN&(!HIF1&(!p38&mTORC1)))))))|((!E2F1&(!ERG&(!ERK&(!FOS&(JUN&(!HIF1&p38))))))|((!E2F1&(!ERG&(!ERK&(FOS&(!HIF1&(!p38&mTORC1))))))|((!E2F1&(!ERG&(!ERK&(FOS&(!HIF1&p38)))))|((!E2F1&(!ERG&(ERK&(!HIF1&(!p38&mTORC1)))))|((!E2F1&(!ERG&(ERK&(!HIF1&p38))))|((!E2F1&(ERG&(!HIF1&(!p38&mTORC1))))|((!E2F1&(ERG&(!HIF1&p38)))|((E2F1&(!ERG&(!ERK&(!FOS&(!JUN&(!GLI&(!HIF1&(!NF_kB&(!SHH&(!TCF&AR_ERG))))))))))|((E2F1&(!ERG&(!ERK&(!FOS&(!JUN&(!GLI&(!HIF1&(!NF_kB&(!SHH&TCF)))))))))|((E2F1&(!ERG&(!ERK&(!FOS&(!JUN&(!GLI&(!HIF1&(!NF_kB&SHH))))))))|((E2F1&(!ERG&(!ERK&(!FOS&(!JUN&(!GLI&(!HIF1&NF_kB)))))))|((E2F1&(!ERG&(!ERK&(!FOS&(!JUN&(GLI&!HIF1))))))|((E2F1&(!ERG&(!ERK&(!FOS&(JUN&!HIF1)))))|((E2F1&(!ERG&(!ERK&(FOS&!HIF1))))|((E2F1&(!ERG&(ERK&!HIF1)))|(E2F1&(ERG&!HIF1))))))))))))))))))))))))))))
        MYC_MAX, (!MXI1&(MYC&(!SMAD&!TGFBR)))
        NCOA3, (p38&!SPOP)
        NCOR1, !AKT
        NCOR2, (!AKT&!FOXO)
        NF1, !PKC
        NF_kB, ((!DNA_Damage&(!E_cadherin&(IKK&(!NCOA3&(!PIP3&Snail)))))|((!DNA_Damage&(!E_cadherin&(IKK&(!NCOA3&PIP3))))|((!DNA_Damage&(!E_cadherin&(IKK&NCOA3)))|((!DNA_Damage&(E_cadherin&(!NCOA3&(!PIP3&Snail))))|((!DNA_Damage&(E_cadherin&(!NCOA3&PIP3)))|((!DNA_Damage&(E_cadherin&NCOA3))|((DNA_Damage&(!NCOA3&(!PIP3&Snail)))|((DNA_Damage&(!NCOA3&PIP3))|(DNA_Damage&NCOA3)))))))))
        NKX3_1, (AR&(!ERG&(PKC&!AR_ERG)))
        Nutrients, Nutrients
        PDK1, ((!HIF1&(!MYC_MAX&PIP3))|((!HIF1&MYC_MAX)|HIF1))
        PHDs, (!Hypoxia&ROS)
        PI3K, ((!EGFR&(!FRS2&TGFBR))|((!EGFR&FRS2)|EGFR))
        PIP3, (!p53&PI3K)
        PKC, (RTK&WNT)
        PTCH1, (GLI&!SHH)
        PTEN, (!NF_kB&(p53&!Snail))
        RAF, ((!EZH2&RAS)|EZH2)
        RAGS, (!Hypoxia&Nutrients)
        RAS, ((!EGFR&(!FRS2&RTK))|((!EGFR&FRS2)|EGFR))
        RB1, (!AMPK&(CHK1_2&(!CyclinB&(!CyclinD&!p53))))
        ROS, ((!COX4I2&(!GSH&Hypoxia))|(COX4I2&!GSH))
        RTK, ((!EGFR&FGFR3)|EGFR)
        RUNX2, ((!ETS1&(!FOXO&(!p38&SMAD)))|((!ETS1&(!FOXO&p38))|(ETS1&!FOXO)))
        Rheb, !TSC1_2
        SHH, !ATR
        SMAD, ((!AR&(!TGFBR&TNFalpha))|(!AR&TGFBR))
        SMO, (!PTCH1&SHH)
        SPOP, SPOP
        Slug, ((!DAXX&(!MDM2&(!NF_kB&(!p53&TCF))))|(!DAXX&(!MDM2&(NF_kB&!p53))))
        Snail, ((!beta_catenin&(!GLI&(!GSK3&(!NF_kB&SMAD))))|((!beta_catenin&(!GLI&(!GSK3&NF_kB)))|(!beta_catenin&(GLI&!GSK3))))
        TAK1, ((!TGFBR&TNFalpha)|TGFBR)
        TCF, ((!beta_catenin&(!ERG&(!TAK1&AR_ERG)))|((!beta_catenin&(ERG&!TAK1))|(beta_catenin&!TAK1)))
        TERT, ((!AKT&(!eEF2&(!GLI&(!HIF1&(!ZBTB17&(!MYC_MAX&(!NF1&(NF_kB&(!p53&!SMAD)))))))))|((!AKT&(!eEF2&(!GLI&(!HIF1&(!ZBTB17&(!MYC_MAX&(NF1&(!p53&!SMAD))))))))|((!AKT&(!eEF2&(!GLI&(!HIF1&(!ZBTB17&(MYC_MAX&(!p53&!SMAD)))))))|((!AKT&(!eEF2&(!GLI&(HIF1&(!ZBTB17&(!p53&!SMAD))))))|((!AKT&(!eEF2&(GLI&(!ZBTB17&(!p53&!SMAD)))))|(AKT&(!eEF2&(!ZBTB17&(!p53&!SMAD)))))))))
        TGFBR, (!MED12&TGFb)
        TGFb, TGFb
        TNFalpha, TNFalpha
        TSC1_2, ((!AKT&(!AMPK&(!ERK&(!HIF1&(!p53&(!p90RSK&!RAF))))))|((!AKT&(!AMPK&(!ERK&(!HIF1&(p53&!p90RSK)))))|((!AKT&(!AMPK&(!ERK&(!HIF1&(p53&(p90RSK&!RAF))))))|((!AKT&(!AMPK&(!ERK&(HIF1&(!p53&!p90RSK)))))|((!AKT&(!AMPK&(!ERK&(HIF1&(!p53&(p90RSK&!RAF))))))|((!AKT&(!AMPK&(!ERK&(HIF1&p53))))|((!AKT&(!AMPK&(ERK&(!HIF1&(p53&(!p90RSK&!RAF))))))|((!AKT&(!AMPK&(ERK&(HIF1&(!p53&(!p90RSK&!RAF))))))|((!AKT&(!AMPK&(ERK&(HIF1&(p53&!p90RSK)))))|((!AKT&(!AMPK&(ERK&(HIF1&(p53&(p90RSK&!RAF))))))|((!AKT&(AMPK&(!ERK&(!HIF1&(!p53&!p90RSK)))))|((!AKT&(AMPK&(!ERK&(!HIF1&(!p53&(p90RSK&!RAF))))))|((!AKT&(AMPK&(!ERK&(!HIF1&p53))))|((!AKT&(AMPK&(!ERK&HIF1)))|((!AKT&(AMPK&(ERK&(!HIF1&(!p53&(!p90RSK&!RAF))))))|((!AKT&(AMPK&(ERK&(!HIF1&(p53&!p90RSK)))))|((!AKT&(AMPK&(ERK&(!HIF1&(p53&(p90RSK&!RAF))))))|((!AKT&(AMPK&(ERK&(HIF1&(!p53&!p90RSK)))))|((!AKT&(AMPK&(ERK&(HIF1&(!p53&(p90RSK&!RAF))))))|((!AKT&(AMPK&(ERK&(HIF1&p53))))|((AKT&(!AMPK&(!ERK&(!HIF1&(p53&(!p90RSK&!RAF))))))|((AKT&(!AMPK&(!ERK&(HIF1&(!p53&(!p90RSK&!RAF))))))|((AKT&(!AMPK&(!ERK&(HIF1&(p53&!p90RSK)))))|((AKT&(!AMPK&(!ERK&(HIF1&(p53&(p90RSK&!RAF))))))|((AKT&(!AMPK&(ERK&(HIF1&(p53&(!p90RSK&!RAF))))))|((AKT&(AMPK&(!ERK&(!HIF1&(!p53&(!p90RSK&!RAF))))))|((AKT&(AMPK&(!ERK&(!HIF1&(p53&!p90RSK)))))|((AKT&(AMPK&(!ERK&(!HIF1&(p53&(p90RSK&!RAF))))))|((AKT&(AMPK&(!ERK&(HIF1&(!p53&!p90RSK)))))|((AKT&(AMPK&(!ERK&(HIF1&(!p53&(p90RSK&!RAF))))))|((AKT&(AMPK&(!ERK&(HIF1&p53))))|((AKT&(AMPK&(ERK&(!HIF1&(p53&(!p90RSK&!RAF))))))|((AKT&(AMPK&(ERK&(HIF1&(!p53&(!p90RSK&!RAF))))))|((AKT&(AMPK&(ERK&(HIF1&(p53&!p90RSK)))))|(AKT&(AMPK&(ERK&(HIF1&(p53&(p90RSK&!RAF))))))))))))))))))))))))))))))))))))))))
        VHL, (!Hypoxia&!ROS)
        WNT, ((!ERG&(!p53&AR_ERG))|(ERG&!p53))
        YWHAZ, AR
        ZBTB17, (AR&!MYC_MAX)
        beta_catenin, ((!AXIN1&(!EZH2&(!GSK3&(!p53&YWHAZ))))|(!AXIN1&(EZH2&(!GSK3&!p53))))
        cFLAR, ((!AKT&(AR&!JNK))|(AKT&!JNK))
        eEF2, !eEF2K
        eEF2K, ((!p70S6kab&p90RSK)|p70S6kab)
        fused_event, fused_event
        mTORC1, (AKT&(!AMPK&!TSC1_2))
        mTORC2, ((!AKT&(!AMPK&(RAGS&Rheb)))|((AKT&(!AMPK&(!RAGS&Rheb)))|((AKT&(!AMPK&RAGS))|(AKT&(AMPK&(RAGS&Rheb))))))
        p14ARF, (E2F1&(MYC_MAX&RAS))
        p15, ZBTB17
        p21, ((!AKT&(!ERK&(!HIF1&(!MDM2&(!ZBTB17&(!MYC_MAX&(!p53&(SMAD&!TERT))))))))|((!AKT&(!ERK&(!HIF1&(!MDM2&(!ZBTB17&(!MYC_MAX&(p53&!TERT)))))))|((!AKT&(!ERK&(!HIF1&(!MDM2&(ZBTB17&(!MYC_MAX&!TERT))))))|(!AKT&(!ERK&(HIF1&(!MDM2&(!MYC_MAX&!TERT))))))))
        p38, (!ERK&(!GADD45&MAP3K1_3))
        p53, ((!Acidosis&(!BCL2&(!CHK1_2&(!HIF1&(!HSPs&(!MDM2&(!p14ARF&p38)))))))|((!Acidosis&(!BCL2&(!CHK1_2&(!HIF1&(!HSPs&(!MDM2&p14ARF))))))|((!Acidosis&(!BCL2&(!CHK1_2&(HIF1&(!HSPs&!MDM2)))))|((!Acidosis&(!BCL2&(CHK1_2&(!HSPs&!MDM2))))|(Acidosis&(!BCL2&(!HSPs&!MDM2)))))))
        p70S6kab, ((!mTORC2&PDK1)|mTORC2)
        p90RSK, ((!ERK&PDK1)|ERK)
    """
    )

    nfvs_mtsNFVS = find_minimum_NFVS(bn_real)

    assert len(nfvs_mtsNFVS) <= 13  # the result of mtsNFVS is 13

    for _i in range(10):
        nfvs = find_minimum_NFVS(bn_real)
        assert nfvs == nfvs_mtsNFVS
