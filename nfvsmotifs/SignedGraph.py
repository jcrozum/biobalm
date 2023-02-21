from __future__ import annotations

import os
import networkx as nx # type: ignore

class SignedGraph:
    V: int
    adjacency_list_positive: dict[str, list[str]]
    adjacency_list_negative: dict[str, list[str]]
    
    def __init__(self, vertex_list: list[str]):
        self.V = len(vertex_list)
        self.adjacency_list_positive = {}
        self.adjacency_list_negative = {}
        for node in vertex_list:
            self.adjacency_list_positive[node] = []
            self.adjacency_list_negative[node] = []
    

    def convert_to_undirected_graph(self) -> nx.DiGraph:
        new_vertex_index = 0

        udGraph = nx.DiGraph()
        for node in self.adjacency_list_positive:
            udGraph.add_node(node)


        for node in self.adjacency_list_negative:
            edgeList = self.adjacency_list_negative[node]

            for v in edgeList:
                udGraph.add_edge(node, v)
                udGraph.add_edge(v, node)

        
        for node in self.adjacency_list_positive:
            edgeList = self.adjacency_list_positive[node]

            for v in edgeList:
                if v != node:
                    new_vertex_index += 1

                    new_vertex = "v_new_" + str(new_vertex_index)
                    udGraph.add_node(new_vertex)
                    udGraph.add_edge(node, new_vertex)
                    udGraph.add_edge(new_vertex, node)
                    udGraph.add_edge(new_vertex, v)
                    udGraph.add_edge(v, new_vertex)


        return udGraph


    def set_edge(self, y: str, x: str, sign: int):
        if sign == 1:
            self.adjacency_list_positive[y].append(x)
        else:
            self.adjacency_list_negative[y].append(x)
    

    def get_self_negative_loops(self) -> list[str]:
        result = []

        for node in self.adjacency_list_negative:
            edgeList = self.adjacency_list_negative[node]

            if node in edgeList:
                result.append(node)

        return result


    def remove_vertex(self, v: str):
        # remove vertex v
        del self.adjacency_list_positive[v]
        del self.adjacency_list_negative[v]
        
        self.V = self.V - 1
        
        # update the edge list of each vertex
        for node in self.adjacency_list_positive:
            if v in self.adjacency_list_positive[node]:
                self.adjacency_list_positive[node].remove(v)

        for node in self.adjacency_list_negative:
            if v in self.adjacency_list_negative[node]:
                self.adjacency_list_negative[node].remove(v)


    def get_negative_degree(self, v: str) -> int:
        neg_deg = 0

        for node in self.adjacency_list_negative:
            if v == node:
                neg_deg += len(self.adjacency_list_negative[node])
            else:
                if v in self.adjacency_list_negative[node]:
                    neg_deg += 1

        return neg_deg
    
    
