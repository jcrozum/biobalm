import os
import networkx as nx # type: ignore

class SignedGraph:
    V: int
    adjacencyListPositive: dict[str, list[str]]
    adjacencyListNegative: dict[str, list[str]]
    
    def __init__(self, listVertex: list[str]):
        self.V = len(listVertex)
        self.adjacencyListPositive = {}
        self.adjacencyListNegative = {}
        for node in listVertex:
            self.adjacencyListPositive[node] = []
            self.adjacencyListNegative[node] = []
    

    def convertToUDGraph(self) -> nx.DiGraph:
        new_vertex_index = 0
        vertexList = []

        udGraph = nx.DiGraph()
        for node in self.adjacencyListPositive.keys():
            vertexList.append(node)
            udGraph.add_node(node)


        for node in self.adjacencyListNegative.keys():
            edgeList = self.adjacencyListNegative[node]

            for v in edgeList:
                udGraph.add_edge(node, v)
                udGraph.add_edge(v, node)

        
        for node in self.adjacencyListPositive.keys():
            edgeList = self.adjacencyListPositive[node]

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


    def setEdge(self, y: str, x: str, sign: int):
        if sign == 1:
            self.adjacencyListPositive[y].append(x)
        else:
            self.adjacencyListNegative[y].append(x)


    def isSelfPositiveLoop(self, v: str) -> bool:
        if not self.adjacencyListPositive.has_key(v):
            print("The vertex does not exist")
            return False

        if v in self.adjacencyListPositive[v]:
            return True
        else:
            return False
    

    def isSelfNegativeLoop(self, v: str) -> bool:
        if not self.adjacencyListNegative.has_key(v):
            print("The vertex does not exist")
            return False

        if v in self.adjacencyListNegative[v]:
            return True
        else:
            return False


    def getSelfPositiveLoops(self) -> list[str]:
        result = []

        for node in self.adjacencyListPositive.keys():
            edgeList = self.adjacencyListPositive[node]

            if node in edgeList:
                result.append(node)

        return result
    

    def getSelfNegativeLoops(self) -> list[str]:
        result = []

        for node in self.adjacencyListNegative.keys():
            edgeList = self.adjacencyListNegative[node]

            if node in edgeList:
                result.append(node)

        return result


    def removeVertex(self, v: str):
        # remove vertex v
        del self.adjacencyListPositive[v]
        del self.adjacencyListNegative[v]
        
        self.V = self.V - 1
        
        # update the edge list of each vertex
        for node in self.adjacencyListPositive.keys():
            if v in self.adjacencyListPositive[node]:
                self.adjacencyListPositive[node].remove(v)

        for node in self.adjacencyListNegative.keys():
            if v in self.adjacencyListNegative[node]:
                self.adjacencyListNegative[node].remove(v)


    def GetDegreeNegative(self, v: str) -> int:
        deg_neg = 0

        for node in self.adjacencyListNegative.keys():
            if v == node:
                deg_neg += len(self.adjacencyListNegative[node])
            else:
                if v in self.adjacencyListNegative[node]:
                    deg_neg += 1

        return deg_neg
    
    
