from __future__ import annotations

import networkx as nx  # type: ignore

from biodivine_aeon import BooleanNetwork # type: ignore
from nfvsmotifs.interaction_graph_utils import find_minimum_NFVS  # type: ignore
from nfvsmotifs.petri_net_translation import network_to_petrinet
from nfvsmotifs.trappist_core import trappist, compute_fixed_point_reduced_STG
from nfvsmotifs.space_utils import percolate_network, percolate_space
from nfvsmotifs.motif_avoidant import detect_motif_avoidant_attractors
from nfvsmotifs.state_utils import state_list_to_bdd


class SuccessionDiagram():
    def __init__(self, network: BooleanNetwork):
        self.network = network
        self.petri_net = network_to_petrinet(self.network)
        nfvs_mtsNFVS = find_minimum_NFVS(self.network)
        retained_set_global = {n: 1 for n in nfvs_mtsNFVS}

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

                # This handles the vacuous case and avoids clingo warnings
                if len(self.G.nodes[sd_node]['fixed_vars']) == self.network.num_vars():
                    maxts = []
                else:
                    maxts = trappist(
                        self.network, problem='max', ensure_subspace=self.G.nodes[sd_node]['fixed_vars'])

                stable_motifs = [x for x in maxts if not (
                    x.items() <= self.G.nodes[sd_node]['fixed_vars'].items())]

                # TODO: merge source motifs

                # TODO: properly create terminal restriction space; for now, just avoid stable motifs
                terminal_restriction_space = ~state_list_to_bdd(stable_motifs)
                nodes = [network.get_variable_name(
                    var) for var in network.variables()]

                retained_set = {}
                for n in retained_set_global:
                    if n not in self.G.nodes[sd_node]['fixed_vars']:
                        retained_set[n] = retained_set_global[n]

                # This handles the vacuous case and avoids clingo warnings
                if len(self.G.nodes[sd_node]['fixed_vars'])+len(retained_set) == self.network.num_vars():
                    retained_set.update(self.G.nodes[sd_node]['fixed_vars'])
                    candidates = [retained_set]
                else:
                    candidates = compute_fixed_point_reduced_STG(
                        self.petri_net, nodes, retained_set,
                        avoid_subspaces=stable_motifs,
                        ensure_subspace=self.G.nodes[sd_node]['fixed_vars'])
                attractors = detect_motif_avoidant_attractors(
                    self.network, self.petri_net,
                    candidates, terminal_restriction_space, AVOIDANCE_ITERATIONS,
                    ensure_subspace=self.G.nodes[sd_node]['fixed_vars'])
                self.G.nodes[sd_node]['attractors'] = attractors

                for fixed_vars in stable_motifs:
                    all_fixed_vars = fixed_vars.copy()
                    all_fixed_vars.update(self.G.nodes[sd_node]['fixed_vars'])
                    fixed_vars_perc, _ = percolate_space(
                        self.network, all_fixed_vars)

                    # TODO: check to see if the reduced network has already been found
                    perc_match_ind = None
                    for sd_node_ind, sd_node_data in self.G.nodes(data=True):
                        if sd_node_data['fixed_vars'] == fixed_vars_perc:
                            perc_match_ind = sd_node_ind
                            break

                    if perc_match_ind is None:
                        ind = self.G.number_of_nodes()
                        self.G.add_node(
                            ind, fixed_vars=fixed_vars_perc.copy(), depth=depth)
                        self.G.add_edge(sd_node, ind, motif=fixed_vars)
                        sd_nodes_at_current_depth.add(ind)
                    else:
                        self.G.add_edge(
                            sd_node, perc_match_ind, motif=fixed_vars)

            sd_nodes_at_last_depth = {
                x for x in sd_nodes_at_current_depth}  # manually copy

        # TODO: compute minimal trap spaces and attractors if we reached max depth or size
