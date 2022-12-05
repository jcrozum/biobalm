'''
This file contains functions used to calculate feedback vertex set for directed graph

The major function is the FVS_local_search(), which calculate maximum sub topological ordering of a graph
The algorithm is given by "Applying local search to feedback vertex set problem"
Author of the paper Philippe Galinier, Eunice Lemamou, Mohamed Wassim Bouzidi
The code mimic the pseudocode given in Page 805 in that paper
The code is written by Gang Yang, Penn State University
'''


import networkx as nx
import random
import math
import array



#the following two function calculate the position for the candidate given the existing topological ordering S
def get_position_minus(candidate_incoming_neighbour, S):
  '''
  get_position_minus return position as just after its numbered in-coming neighbours
  As in the paper, the function return i_minus(v)
  '''
  position = 1
  for x in range(len(S)-1,-1,-1):             #we loop through index from the end as we try to find the largest index of numbered incoming neighbour
    if S[x] in candidate_incoming_neighbour:
      position = x+2                          #2 comes from we need to put candidate after the incoming neighbour and also the position count from 1 instead of 0
      return position
  return position                             #put the candidate in the first position if there is no numbered incoming neighbour



def get_position_plus(candidate_outgoing_neighbour, S):
  '''
  get_position_plus return position as just before its numbered out-going neighbours
  As in the paper, the function return i_plus(v)
  '''
  position = 1+len(S)
  for x in range(len(S)):                     #we loop through index from the beginning as we try to find the smallest index of numbered outgoing neighbour
    if S[x] in candidate_outgoing_neighbour:
      position = x+1                          #1 comes from the fact position count from 1 instead of 0
      return position
  return position                             #put the candidate in the first position if there is no numbered outgoing neighbour



def FVS_local_search(G_input, T_0, alpha, maxMvt, maxFail, randomseed=None):
  '''
  Returns an maximum sub topological ordering of a DiGraph G.
  FVS is G_input \ the topological ordering
  A topological ordering of a graph is an ordering of vertices such that
  the starting point of every arc occurs earlier in the ordering than
  the endpoint of the art.

  This algorithm is a fast approximation based on simulating annealing(SA) of a noval local search strategy [1]_.
  The code follows the pseudocode given in Page 805 in that paper.


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
  An approximation of the maximum ordering of the given graph as a list.


  Notes
  -----
  The code is written by Gang Yang, Department of Physics, Penn State University


  References
  ----------
  ..[1] Galinier, P., Lemamou, E. & Bouzidi, M.W. J Heuristics (2013) 19: 797. doi:10.1007/s10732-013-9224-z

  '''
  #set a random seed
  random.seed(randomseed)
  G=G_input.copy()
  N= len(G.nodes())           #number of nodes in the graph
  edges=list(G.edges())


  #Initialization
  T = T_0                     #set initial temperatue
  nbFail = 0                  #Outer loop counter to record failure times
  S = []                      #list to record the ascending ordering
  S_optimal = []              #list to record the optimal ordering


  #calculate parent and child node for each node
  parent = [{} for i in range(N)]
  child = [{} for i in range(N)]

  for i in range(len(edges)):
    (a, b) = edges[i]
    child[int(a)][int(b)]=None
    parent[int(b)][int(a)]=None

  #all the nodes that are self_loops
  self_loops = [edges[i][0] for i in range(len(edges)) if edges[i][0]==edges[i][1]]
  #all the nodes that is not in S
  unnumbered=[x for x in range(N) if x not in self_loops]
  N_unnumbered = len(unnumbered)
  while nbFail< maxFail:
    nbMvt = 0       #Inner loop counter to record movement times
    failure = True  #if cardinal of S increase after one inner loop, failure will be set to false
    while nbMvt < maxMvt:
      candidate_index = random.randint(0, N_unnumbered-1)
      candidate = unnumbered[candidate_index]         #random pick a node from all unnumbered node
      position_type = random.randint(0,1)             #random choose a position type
      #calculate incoming/outgoing neighbour for the candidate node, store as keys in the dict
      candidate_incoming_neighbour = parent[candidate]
      candidate_outgoing_neighbour = child[candidate]
      #see how to find the position on Page 803 of the paper
      #position_type=1 means just after incoming neighbours
      if position_type==1:
        position = get_position_minus(candidate_incoming_neighbour,S)
      #position_type=0 means just before outgoint neighbours
      elif position_type==0:
        position = get_position_plus(candidate_outgoing_neighbour,S)

      #first, insert the candidate to the given position
      S_trail = S[:]  #copy the list
      S_trail.insert(position-1,candidate)

      #second remove all the conflict
      #break the sequence into two parts: before and after the candidate node and
      S_trail_head= S_trail[:position-1]
      S_trail_tail= S_trail[position:]
      #determine conflict node See page 801
      if position_type==1:
        CV_pos=[]    #conflict before the newly inserted node in the topological ordering
        for x in range(len(S_trail_head)):
          nodetemp=S_trail_head[x]
          #print nodetemp,candidate_outgoing_neighbour,nodetemp in candidate_outgoing_neighbour
          if nodetemp in candidate_outgoing_neighbour:
            CV_pos.append(nodetemp)
        conflict=CV_pos   #there won't be conflict after the inserted node as the node inserted after its incoming neighbour
      elif position_type==0:
        CV_neg=[]    #conflict after the newly inserted node in the topological ordering
        for x in range(len(S_trail_tail)):
          nodetemp=S_trail_tail[x]
          #print nodetemp,candidate_incoming_neighbour,nodetemp in candidate_incoming_neighbour
          if nodetemp in candidate_incoming_neighbour:
            CV_neg.append(nodetemp)
        conflict=CV_neg   #there won't be conflict before the inserted node as the node inserted before its outgoing neighbour
      #finally remove the conflict node
      N_conflict=len(conflict)
      if N_conflict>0:
        for i in range(N_conflict):
          S_trail.remove(conflict[i])

      #third, evaluate the move
      #delta_move=-len(S_trail)+len(S)
      delta_move=N_conflict-1
      #accept all the move that is not harmful, otherwise use metrospolis algorithm
      if delta_move<=0 or math.exp(-delta_move/float(T))>random.random():
        S = S_trail[:]
        #update unnumbered nodes
        unnumbered.remove(candidate)   #remove the node just inserted
        if N_conflict>0:               #add all the conflict node just removed
          for i in range(N_conflict):
            unnumbered.append(conflict[i])
        N_unnumbered+=delta_move
        nbMvt = nbMvt+1
        #update S_optimal only when there is increase in cardinal of the sequence
        if len(S)>len(S_optimal):
          S_optimal = S[:]
          failure = False
        if N_unnumbered==0:
          return S_optimal

    #Increment the failure times if no progress in size of the sequence
    if failure==True:
      nbFail+=1
    else:    #otherwise reset the num of failure times
      nbFail=0
    #shrink the temperatue by factor alpha
    T=T*alpha
    #print T
    #print nbFail
  return S_optimal



