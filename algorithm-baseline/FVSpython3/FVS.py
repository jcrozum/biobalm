'''
This file contains functions used to calculate feedback vertex set for graph/directed graph.
'''

import networkx as nx
#import FVS_localsearch_10_cython as FVS10
import FVSpython3.FVS_localsearch_10_python as FVS10
#import FVS_localsearch_10_python as FVS10

def python_format(Ginput):
  '''
  The FVS_local_search function are written based on the fact that the node are named as
  sequential integer number from 0 to N-1. This function convert the input graph into a
  graph that follows this criteria and can be used by FVS_local_search. The function returns
  the formatted graph, this mapping from original name to new name and its inverse mapping.
  '''
  G=nx.DiGraph()
  N=len(G.nodes())
  i=0
  mapping={}
  inverse_mapping={}
  for node in Ginput.nodes():
    mapping[node]=i
    inverse_mapping[i]=node
    i+=1
  new_edgelist=[(mapping[edge[0]],mapping[edge[1]]) for edge in Ginput.edges()]
  G.add_edges_from(new_edgelist)
  return G,mapping,inverse_mapping


def FVS(G, T_0 = 0.6, alpha = 0.99, maxMvt_factor = 5, maxFail = 50, randomseed=None):
  '''
  Returns an approximation of minimum Feedback Vertex Set (FVS) for a DiGraph G.
  An undirected graph will be converted to a directed graph.
  If the input graph is weakly connected, FVS_weakly_connected() will be called.
  If the input graph is not weakly connected, i.e. have disconnected components,
  FVS_disconnected() will be called.


  From Wikipedia:
  A Feedback Vertex Set of a graph is a set of vertices whose removal leaves a graph without cycles.
  In other words, each Feedback Vertex Set contains at least one vertex of any cycle in the graph.
  The Feedback Vertex Set problem for directed graph is an NP-complete problem in computational complexity theory.

  This algorithm is a fast approximation based on simulating annealing(SA) of a noval local search strategy [1]_.
  The function mainly called a function named FVS_local_search.
  FVS_local_search calcualte the longest sub topological ordering of G_input (directed graph)
  based on simulated annealing(SA) of local search of topological ordering.
  FVS is then G \ the topological ordering.
  The code follows the pseudocode given in Page 805 in that paper.
  This algorithm is claimed to outperform the existing heuristic algorithm,
  greedy adaptive search procedure (GRASP) by Pardalos et al.

  Parameters
  ----------
  G : NetworkX Graph/DiGraph, result for MultiGraph is not tested
  T_0 : the initial temperature  in SA
  alpha : the cooling multiplier for temperatue in the geometric cooling regime
  maxMvt_factor : maxMvt_factor times network size is the number of iterations for the inner loop given a fixed temperatue
  maxFail : FVS_local_search stops when maxFail number of outloops (temperatue) fail to improve the result
  randomseed: random seed for generating random numbers

  Returns
  -------
  An approximation of the minimum FVS of the given graph as a list.


  Examples
  --------
  #calculate a FVS for a random graph
  >>>G1 = nx.gnm_random_graph(6, 12, seed=None, directed=True)
  >>>G1_FVS = FVS.FVS(G1)
  #calculate a FVS for a toy example
  >>>G2 = nx.DiGraph()
  >>>G2.add_edges_from([('A','B'),('B','C'),('C','A'),('A','D'),('D','A')])
  >>>G2_FVS = FVS.FVS(G2)
  #result should be ['A']

  Notes
  -----
  The code is written by Gang Yang, Department of Physics, Penn State University

  References
  ----------
  ..[1] Galinier, P., Lemamou, E. & Bouzidi, M.W. J Heuristics (2013) 19: 797. doi:10.1007/s10732-013-9224-z

  '''
  #make sure input is a directed graph, undirected graph will be converted into a directed graph
  if not nx.is_directed(G):
    print("Warning: undirected graph is converted to directed graph!")
    Gtemp=G.to_directed()
  else:
    Gtemp=G.copy()
  if nx.is_weakly_connected(Gtemp):
    result = FVS_weakly_connected(Gtemp,T_0,alpha,maxMvt_factor,maxFail,randomseed)
  else:
    result = FVS_disconnected(Gtemp,T_0,alpha,maxMvt_factor,maxFail,randomseed)
  return result




def FVS_weakly_connected(G, T_0, alpha, maxMvt_factor, maxFail,randomseed):
  '''
  Returns an approximation of minimum Feedback Vertex Set (FVS) for a weakly connected DiGraph G.
  This function is part of FVS(). Do not call this function directly.
  See more information at FVS().

  Parameters
  ----------
  G : NetworkX DiGraph (required to be weakly connected)
  T_0 : the initial temperature  in SA
  alpha : the cooling multiplier for temperatue in the geometric cooling regime
  maxMvt_factor : maxMvt_factor times network size is the number of iterations for the inner loop given a fixed temperatue
  maxFail : FVS_local_search stops when maxFail number of outloops (temperatue) fail to improve the result
  randomseed: random seed for generating random numbers

  Returns
  -------
  An approximation of the minimum FVS of the given graph as a list.


  Notes
  -----
  The code is written by Gang Yang, Department of Physics, Penn State University

  '''
  #test whether it is already a Directed Acyclic Graph
  self_loops = [edge[0] for edge in G.edges() if edge[0]==edge[1] ]   #all the nodes that are self_loops
  Gtemp = G.copy()
  Gtemp.remove_nodes_from(self_loops)
  if nx.is_directed_acyclic_graph(Gtemp):
    return self_loops
  #map the network into a format required by FVS_local search
  G_formatted, mapping, inverse_mapping = python_format(G)
  N = len(G.nodes())
  #maximum iteration times for inner loop
  maxMvt = maxMvt_factor*N
  #calculate the longest topological order inside the given graph
  S_optimal_formatted = FVS10.FVS_local_search(G_formatted, T_0, alpha, maxMvt, maxFail, randomseed)
  #map the topological order back to the names in the original graph
  S_optimal = [inverse_mapping[node] for node in S_optimal_formatted]
  #minimum FVS is all the nodes minus nodes in the topological order
  FVS_set = set(G.nodes())-set(S_optimal)
  FVS = [node for node in FVS_set]
  #check that we do obtain a FVS through checking G\FVS is a Directed Acyclic Graph
  Gprime=G.copy()
  Gprime.remove_nodes_from(FVS)
  assert nx.is_directed_acyclic_graph(Gprime)
  return FVS


def FVS_disconnected(G, T_0, alpha, maxMvt_factor, maxFail, randomseed):
  '''
  Returns an approximation of minimum Feedback Vertex Set (FVS) for a disconnected DiGraph G.
  This function is part of FVS(). Do not call this function directly.
  See more information at FVS().

  For each weakly connnected component of the original graph,
  call FVS_weakly_connected to calculate FVS of each component.
  The result is the union of the FVS of each component.


  Parameters
  ----------
  G : A disconnected NetworkX DiGraph
  T_0 : the initial temperature  in SA
  alpha : the cooling multiplier for temperatue in the geometric cooling regime
  maxMvt_factor : maxMvt_factor times network size is the number of iterations for the inner loop given a fixed temperatue
  maxFail : FVS_local_search stops when maxFail number of outloops (temperatue) fail to improve the result
  randomseed: random seed for generating random numbers

  Returns
  -------
  An approximation of the minimum FVS of the given graph as a list.


  Notes
  -----
  The code is written by Gang Yang, Department of Physics, Penn State University
  '''

  wccs = (G.subgraph(c).copy() for c in nx.weakly_connected_components(G))
  #wccs = nx.weakly_connected_component_subgraphs(G)
  FVS=[]
  for wcc in wccs:
    FVS_temp = FVS_weakly_connected(wcc, T_0, alpha, maxMvt_factor, maxFail, randomseed)
    FVS.extend(FVS_temp)
  return FVS


