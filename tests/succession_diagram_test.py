from nfvsmotifs.SuccessionDiagram import SuccessionDiagram
from biodivine_aeon import BooleanNetwork, SymbolicAsyncGraph, find_attractors # type: ignore
import sys

def test_succession_diagram_structure():
    bn = BooleanNetwork.from_bnet("""
        x1, x2
        x2, x1
        x3, !x3
        """)

    SD = SuccessionDiagram(bn)
    SD.expand_node(SD.root(), depth_limit=None)
    assert SD.G.number_of_nodes() == 3
    assert SD.G.number_of_edges() == 2
    assert max(d['depth'] for n,d in SD.G.nodes(data=True)) == 1
    
    
# TODO: add tests for a wider variety of networks

def test_attractor_detection(network_file):
    # TODO: Once attractor detection is faster, we should increase this limit.
    # Right now, checking attractors in larger succession diagrams would often time out our CI.
    NODE_LIMIT = 100

    # This is unfortunately necessary for PyEDA Boolean expression parser (for now).
    sys.setrecursionlimit(150000)

    # TODO: Remove these once method is fast enough.
    print(network_file)
    if network_file.endswith("146.bnet"):
        # For this model, we can compute the 100 SD nodes, but it takes a very long time
        # and the SD is larger, so we wouldn't get to attractor computation anyway.
        NODE_LIMIT = 10

    bn = BooleanNetwork.from_file(network_file)
    bn = bn.infer_regulatory_graph()
    stg = SymbolicAsyncGraph(bn)

    # Compute the succession diagram.
    sd = SuccessionDiagram(bn)
    expanded = sd.expand_node(sd.root(), depth_limit=1000, node_limit=NODE_LIMIT)

    # SD must be fully expanded, otherwise we may miss some results.
    # If SD is not fully expanded, we just skip this network.
    if expanded >= NODE_LIMIT:
        return
    
    print(expanded)

    # TODO: Remove these once method is fast enough. 
    if network_file.endswith("075.bnet"):
        # It seems that with current NFVS, the clingo fixed-point part takes too long. There are
        # better NFVS-es that we could try, but we first need to make the NFVS algorithm deterministic.
        return


    # Compute attractors in diagram nodes.
    # TODO: There will probably be a method that does this in one "go".
    nfvs_attractors = []
    for i in range(sd.G.number_of_nodes()):
        attr = sd.expand_attractors(i)
        for a in attr:
            # Just a simple sanity check.
            assert len(a) == bn.num_vars()
        if len(attr) > 0:
            nfvs_attractors += attr

    # Compute symbolic attractors using AEON.
    symbolic_attractors = find_attractors(stg)

    # Check that every "seed" returned by SuccessionDiagram appears in
    # some symbolic attractor, and that every symbolic attractor contains
    # at most one such "seed" state.
    for seed in nfvs_attractors:
        symbolic_seed = stg.fix_subspace({ k:bool(v) for k,v in seed.items() })
        found = None
        
        # The "seed" state must have a symbolic attractor (and that
        # attractor mustn't have been removed yet).
        for i in range(len(symbolic_attractors)):
            if symbolic_seed.is_subset(symbolic_attractors[i]):
                found = i
        assert found is not None

        symbolic_attractors.pop(found)

    print("Attractors:", len(nfvs_attractors))

    # All symbolic attractors must be covered by some seed at this point.
    assert len(symbolic_attractors) == 0
