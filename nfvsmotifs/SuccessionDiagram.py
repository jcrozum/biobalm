from __future__ import annotations

import networkx as nx # type: ignore

from biodivine_aeon import BooleanNetwork  # type: ignore
from nfvsmotifs.petri_net_translation import network_to_petrinet
from nfvsmotifs.trappist_core import trappist, compute_fixed_point_reduced_STG
from nfvsmotifs.space_utils import percolate_network, percolate_space
from nfvsmotifs.motif_avoidant import detect_motif_avoidant_attractors
from nfvsmotifs.state_utils import state_list_to_bdd

class SuccessionDiagram():
    def __init__(self, network: BooleanNetwork):
        self.network = network

        AVOIDANCE_ITERATIONS = 10
        MAX_DEPTH = float('inf')
        MAX_SIZE = float('inf')
        depth = 0

        self.G = nx.DiGraph()
        self.G.add_node(0, fixed_vars={}, depth=depth)
        sd_nodes_at_last_depth = {0}

        while (sd_nodes_at_last_depth and depth < MAX_DEPTH and self.G.number_of_nodes() < MAX_SIZE):
            depth += 1
            sd_nodes_at_current_depth = set()
            for sd_node in sd_nodes_at_last_depth:
                reduced_network = percolate_network(
                    self.network, self.G.nodes[sd_node]['fixed_vars'])

                stable_motifs = [x for x in trappist(reduced_network, problem='max') if not(
                    x.items() <= self.G.nodes[sd_node]['fixed_vars'].items())]
                
                # TODO: merge source motifs

                # TODO: properly create terminal restriction space; for now, just avoid stable motifs
                terminal_restriction_space = ~state_list_to_bdd(stable_motifs)
                petri_net = network_to_petrinet(reduced_network)
                nodes = [network.get_variable_name(var) for var in network.variables()]
                retained_set = {n:1 for n in nodes}
                
                # TODO: properly compute candidates
                candidates = compute_fixed_point_reduced_STG(petri_net, nodes, retained_set, avoid_subspaces = stable_motifs)
                attractors = detect_motif_avoidant_attractors(
                    reduced_network, petri_net, candidates, terminal_restriction_space, AVOIDANCE_ITERATIONS)

                self.G.nodes[sd_node]['attractors'] = attractors

                for fixed_vars in stable_motifs:
                    fixed_vars_perc, _ = percolate_space(
                        reduced_network, fixed_vars)
                    
                    # TODO: check to see if the reduced network has already been found
                    perc_match_ind = None

                    if perc_match_ind is None:
                        ind = self.G.number_of_nodes()
                        self.G.add_node(
                            ind, fixed_vars=fixed_vars_perc, depth=depth)
                        self.G.add_edge(sd_node, ind, motif=fixed_vars)
                        sd_nodes_at_current_depth.add(ind)
                    else:
                        self.G.add_edge(
                            sd_node, perc_match_ind, motif=fixed_vars)

            sd_nodes_at_last_depth = {
                x for x in sd_nodes_at_current_depth}  # manually copy

        # TODO: compute minimal trap spaces and attractors if we reached max depth or size
