# README

# I)	THE METHOD
This repository contains functions used to approximate the minimum feedback vertex set (FVS) of a directed graph, a problem whose exact solution is known to be NP-hard. This algorithms was used to estimate the size of the minimum FVS in a recent paper "Zañudo, J. G. T., Yang, G., & Albert, R. (2017). PNAS 114: 28, 7234–7239, doi: 10.1073/pnas.1617387114".

This algorithm identifies a near-minimum feedback vertex set using simulated annealing (SA) and a local search of topological ordering. The algorithm is describe in the paper "Galinier, P., Lemamou, E. & Bouzidi, M.W. J Heuristics (2013) 19: 797. doi:10.1007/s10732-013-9224-z". The code follows the pseudocode given in Page 805 in that paper.

The originally code was written in Python 2.7, and this updated version seems to work for Python 3 (although it is still not fully tested). The module requires NetworkX 1.11.

Another version written in Cython is available at https://github.com/yanggangthu/FVS_cython. The cython version is typically 5-10 times faster than the python version. It has not been updated for use in Python 3.


# II) STRUCTURE OF MODULE
Related functions are stored in FVS.py and FVS_localsearch_10_python.py. 

Core functions of finding a maximum sub topological ordering of a given graph (which is equivalent to identifying a minimum FVS) is written in FVS_localsearch_10_python.py.

FVS.py is a wrapper, that deals with graph input and identifies a near-minimum FVS by subtracting the nodes in the topological ordering from all the nodes.

FVS_test.py contains three examples illustrating how to use the code.

# III)	INSTRUCTIONS  

To use this module, import FVS and call FVS() just as a regular function. The function can take 6 paramters listed below and only the first (the graph) is neccessary.

This was tested in anaconda3 creating an environment "py37_FVS" with Python=3.7 and networkx=1.1:

```
conda create -n py3_FVS python=3.7 networkx=1.11
```

From there, the FVS_python3 folder with the contents of the repository was moved to *anaconda3/envs/py37_FVS/lib/python3.7/site-packages*. With that, the *FVS_test.py* runs as it should:

```
python FVS_test.py

[(0, 4), (0, 1), (1, 0), (2, 1), (2, 4), (2, 5), (2, 3), (3, 5), (3, 4), (4, 1), (5, 2), (5, 0)]
[0, 2]
[('A', 'B'), ('A', 'D'), ('B', 'C'), ('C', 'A'), ('D', 'A')]
['A']
[('A', 'B'), ('B', 'C'), ('C', 'A')]
['B']
```

Parameters
----------
G : NetworkX Graph/DiGraph, result for MultiGraph is not tested  
T_0 : the initial temperature in SA  
alpha : the cooling multiplier for temperatue in the geometric cooling regime  
maxMvt_factor : maxMvt_factor times network size is the number of iterations for the inner loop given a fixed temperatue  
maxFail : FVS_local_search stops when maxFail number of outloops (temperatue) fail to improve the result  
randomseed: random seed for generating random numbers  

Default Parameter Values  
-----------------------
T_0 = 0.6, alpha = 0.99, maxMvt_factor = 5, maxFail = 50, randomseed=None  
The default values are suggested by the author of the paper.  
T_0 and maxFail are chosen after a limited number of preliminary experiments
alpha is chosen more arbitrarily, however alpha should be a positive number slightly smaller than 1.  
Increase alpha or maxMvt_factor or maxFail will increase the time of finding FVS.

Returns
-------
An approximation of the minimum FVS of the given graph as a list.

# IV) EXAMPLES

```
import networkx as nx  
from FVS_python3 import FVS as FVS
```
Here we construct an example with an optimal solution. G2_FVS shoule be ['A'] as a list.

```  
G2=nx.DiGraph()  
G2.add_edges_from([('A','B'),('B','C'),('C','A'),('A','D'),('D','A')])  
G2_FVS=FVS.FVS(G2)
```

Here we construct an example of three-node feedback loops.   
We show how you change all the parameters and set a random seed.  
Your result should be the same with the same randomseed.  
```  
G3=nx.DiGraph()  
G3.add_edges_from([('A','B'),('B','C'),('C','A')])   
G3_FVS=FVS.FVS(G3, T_0=0.6, alpha=0.99, maxMvt_factor=5, maxFail=50, randomseed=1) 
``` 



# V)	COPYRIGHT


The MIT License (MIT)

Copyright (c) 2017 Gang Yang and Reka Albert.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
