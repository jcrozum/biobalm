#### The case of source nodes

Consider the BN

    A, A
    B, B
    C, A & B
    D, D | A

Max. trap spaces at level 0 by running `for x in ar.succession_diagram.motif_reduction_dict[0].stable_motifs: print(x)`

    {'D': 1}
    {'B': 1}
    {'B': 0}
    {'A': 1}
    {'A': 0}

Get the succession diagram by running

    import pystablemotifs as sm
    import pyboolnet
    import pystablemotifs.export as ex
    import networkx as nx
    max_simulate_size=20
    primes = sm.format.import_primes('test.bnet')
    ar = sm.AttractorRepertoire.from_primes(primes, max_simulate_size=max_simulate_size)
    GM=ex.networkx_succession_diagram_motif_based(ar,include_attractors_in_diagram=True)
    ex.plot_nx_succession_diagram(GM)

However, in the succession diagram, the roots are `{A = 0, B = 0}`, `{A = 0, B = 1}`, `{A = 1, B = 0}`, and `{A = 1, B = 1}`.
Why?

One possible improvement is to add constraints to the encoded ASP to ensure that **in every max. trap space, a source node is always a fixed node**.
This still preserves the correctness of nfvs-motifs (e.g., `{'D': 1}` still appears later in the succession diagram).
This may also reduce the running time of nfvs-motifs.