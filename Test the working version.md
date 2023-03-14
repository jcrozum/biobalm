# nfvs-motifs: Test the working version

## An infinite loop in Succession Diagram construction

I see that this error has been fixed by Jordan.

Jordan uses 
		
		retained_set = self.G.nodes[sd_node]['fixed_vars']

It seems that you skip the motif-avoidant attractor checking because the candidate set will have no state.

## Some changes

### Properly set the retained_set

	nfvs_mtsNFVS = find_minimum_NFVS(self.network)
	retained_set_global = {n:1 for n in nfvs_mtsNFVS}
	...
	retained_set = {}
	for n in retained_set_global:
	    if n not in self.G.nodes[sd_node]['fixed_vars']:
	       retained_set[n] = retained_set_global[n]

Now, we can properly compute the candidate set, but this set may contain more states because `b_i` is always set to 1.

### Properly compute terminal_restriction_space

Regarding avoiding stable motifs, the current code seems to be incorrect because it does not consider the covered sub-space.

	terminal_restriction_space = ~state_list_to_bdd(stable_motifs)

It should be

	terminal_restriction_space = ~state_list_to_bdd(stable_motifs) & state_list_to_bdd([self.G.nodes[sd_node]['fixed_vars']])

### Optimize the function detect_motif_avoidant_attractors

	def detect_motif_avoidant_attractors(
	    network: BooleanNetwork,
	    petri_net: DiGraph,
	    candidates: list[dict[str, int]],
	    terminal_restriction_space: BinaryDecisionDiagram,
	    max_iterations: int
	) -> list[dict[str, int]]:
	    if len(candidates) == 0:
	        return []
	
	    candidates = _preprocess_candidates(network, candidates, terminal_restriction_space, max_iterations)
	
	    if len(candidates) == 0:
	        return []
	        
	    return _filter_candidates(petri_net, candidates, terminal_restriction_space)

The current code is common for both the following cases.
The first one is to detect motif-avoidant attractors.
The second one is to detect attractors inside minimal trap spaces.
For the first case, the best case of `_preprocess_candidates` is to exclude all states from `candidates`, i.e., `len(candidates) = 0` after finishing `_preprocess_candidates`.
For the second case, the best case of `_preprocess_candidates` is that `len(candidates) = 1` after finishing `_preprocess_candidates`.
I add the `is_in_an_mts` to separate the above cases.

	def detect_motif_avoidant_attractors(
		    network: BooleanNetwork,
		    petri_net: DiGraph,
		    candidates: list[dict[str, int]],
		    terminal_restriction_space: BinaryDecisionDiagram,
		    max_iterations: int,
		    is_in_an_mts: bool
		) -> list[dict[str, int]]:
		    if len(candidates) == 0:
		        return []
	
		    if len(candidates) == 1 and is_in_an_mts:
		        return candidates
		    
		    candidates = _preprocess_candidates(network, candidates, terminal_restriction_space, max_iterations)
		    
		    if len(candidates) == 0:
		        return []
	
		    if len(candidates) == 1 and is_in_an_mts:
		        return candidates
	
		    return _filter_candidates(petri_net, candidates, terminal_restriction_space)

In the file `SuccessionDiagram.py`, we call

	attractors = detect_motif_avoidant_attractors(
	                reduced_network, petri_net, candidates, 
	                terminal_restriction_space, AVOIDANCE_ITERATIONS, 
	                is_in_an_mts = (len(stable_motifs) == 0)
	             )


## Error with Pint-reach

	def detect_motif_avoidant_attractors(
			    network: BooleanNetwork,
			    petri_net: DiGraph,
			    candidates: list[dict[str, int]],
			    terminal_restriction_space: BinaryDecisionDiagram,
			    max_iterations: int
			) -> list[dict[str, int]]:
			    if len(candidates) == 0:
			        return []
		    
			    candidates = _preprocess_candidates(network, candidates, terminal_restriction_space, max_iterations)
	
			    if len(candidates) == 0:
			        return []
			        
			    return _filter_candidates(petri_net, candidates, terminal_restriction_space)

If using the current code of the function `detect_motif_avoidant_attractors` and the new code for properly computing `terminal_restriction_space`, the function `_filter_candidates` is invoked.
But I get an error when `_Pint_reachability` is invoked.
The example Boolean network is

	x1, x2
	x2, x1
	x3, !x3

And the error is

	...
	File "/home/giang-trinh/git_repo/nfvs-motifs/nfvsmotifs/motif_avoidant.py", line 178, in _Pint_reachability
	    pint_model = InMemoryModel(petri_net_as_automata_network(petri_net))
	  File "/home/giang-trinh/anaconda3/lib/python3.9/site-packages/pypint/model.py", line 311, in __init__
	    self.load()
	  File "/home/giang-trinh/anaconda3/lib/python3.9/site-packages/pypint/model.py", line 227, in load
	    raise e
	pypint.tools.PintProcessError: Command 'pint-export -l nbjson' returned non-zero exit status 2
	Fatal error: exception Failure("Line 68 char 6: Syntax error")

Samuel: Could you please look at this error?

## Test some real-world models

`pystablemotifs`: `max_sim = 20`
`nfvs-motifs`: the updated version

| **Network** | **n** | **# SD nodes** | **pystablemotifs (s)** | **nfvs-motifs (s)** |
| -------------- | ----: | ----: | --------: | -------: | --------: |
| T-LGL-SURVIVAL-NETWORK-2008 | 58 | 2737 | > 600 | 202.43 |
| MELANOGENESIS | 52 | 3 | N/A (1.43) | 33.48 |
| COLITIS-ASSOCIATED-COLON-CANCER | 66 | 85 | 1895.23 | 5.92 |

`# SD nodes` is obtained from the result of `nfvs-motifs`.
This result is not optimized because we have not implemented the feed-forward detection yet.
`N/A` means that the attractor existence has not been fully explored.

For the `MELANOGENESIS` network, `pystablemotifs` is faster than `nfvs-motifs`.
The reason is that `nfvs-motifs` spends most of its running time to find attractors inside minimal trap spaces, which is not performed with `pystablemotifs`.
Moreover, `nfvs-motifs` has not been optimized yet.

I think these initial results are promising.
We need to test and optimize the current code more.





***

## {2023-03-14}Test the working version

### PyBoolNet's repository

I set `AVOIDANCE_ITERATIONS = 20000`.
Excluded two models with too many source nodes.

	arellano_rootstem.bnet (9 nodes)
	# attractors = 4
	# nodes = 7, time = 0.03
	=========
	calzone_cellfate.bnet (28 nodes)
	# attractors = 27
	# nodes = 53, time = 0.96
	=========
	dahlhaus_neuroplastoma.bnet (23 nodes)
	# attractors = 32
	# nodes = 49, time = 200.25
	=========
	davidich_yeast.bnet (10 nodes)
	# attractors = 12
	# nodes = 47, time = 0.41
	=========
	dinwoodie_life.bnet (15 nodes)
	# attractors = 7
	# nodes = 37, time = 0.52
	=========
	dinwoodie_stomatal.bnet (13 nodes)
	# attractors = 1
	# nodes = 2, time = 0.05
	=========
	faure_cellcycle.bnet (10 nodes)
	# attractors = 2
	# nodes = 4, time = 2.82
	=========
	grieco_mapk.bnet (53 nodes)
	# attractors = 18
	# nodes = 43, time = 458.60
	=========
	irons_yeast.bnet (18 nodes)
	# attractors = 1
	# nodes = 1, time = 15.86
	=========
	klamt_tcr.bnet (40 nodes)
	# attractors = 8
	# nodes = 9, time = 0.21
	=========
	krumsiek_myeloid.bnet (11 nodes)
	# attractors = 6
	# nodes = 18, time = 0.12
	=========
	multivalued.bnet (13 nodes)
	# attractors = 4
	# nodes = 7, time = 0.02
	=========
	n12c5.bnet (12 nodes)
	# attractors = 5
	# nodes = 15, time = 0.24
	=========
	n3s1c1a.bnet (3 nodes)
	# attractors = 2
	# nodes = 3, time = 0.00
	=========
	n3s1c1b.bnet (3 nodes)
	# attractors = 2
	# nodes = 3, time = 0.00
	=========
	n5s3.bnet (5 nodes)
	# attractors = 3
	# nodes = 5, time = 0.03
	=========
	n6s1c2.bnet (6 nodes)
	# attractors = 3
	# nodes = 6, time = 0.02
	=========
	n7s3.bnet (7 nodes)
	# attractors = 3
	# nodes = 6, time = 0.02
	=========
	raf.bnet (3 nodes)
	# attractors = 2
	# nodes = 3, time = 0.00
	=========
	randomnet_n15k3.bnet (15 nodes)
	# attractors = 3
	# nodes = 5, time = 0.09
	=========
	randomnet_n7k3.bnet (7 nodes)
	# attractors = 10
	# nodes = 27, time = 0.17
	=========
	remy_tumorigenesis.bnet (35 nodes)
	# attractors = 25
	# nodes = 63, time = 153.32
	=========
	saadatpour_guardcell.bnet (13 nodes)
	# attractors = 1
	# nodes = 2, time = 0.05
	=========
	tournier_apoptosis.bnet (12 nodes)
	# attractors = 3
	# nodes = 6, time = 0.07
	=========
	xiao_wnt5a.bnet (7 nodes)
	# attractors = 4
	# nodes = 7, time = 0.02
	=========
	zhang_tlgl.bnet (60 nodes)
	# attractors = 156
	# nodes = 418, time = 2024.26
	=========
	zhang_tlgl_v2.bnet (60 nodes)
	_filter_candidates
	# attractors = 258
	# nodes = 957, time = 3055.04
	=========

Observations:
+ First, I think all the issues I reported have been addressed.
+ Second, The number of attractors returned is correct for every model.
+ Third, however, the running time is still **large** for some models (e.g., `dahlhaus_neuroplastoma`, `grieco_mapk`, `remy_tumorigenesis`, `zhang_tlgl`, and `zhang_tlgl_v2`).
+ For the `zhang_tlgl_v2` model, it needs to call the `_filter_candidates` function. From my experience with `mtsNFVS`, we can avoid calling the `_filter_candidates` function for most real-world models (if Preprocessing SSF is good enough).

The reason for the slow running time may be the fact that `nfvs-motifs` has been not optimized.
Hence, I think our next step is to optimize `nfvs-motifs`.



### BBM repository

Sam has also tested the current code on the BBM repository.

For every but one model where the SD is successfully obtained, the number of attractors is the same as that obtained by AEON.

Need to investigate this model more.



## Optimization

### Improve Preprocessing SSF

	def _preprocess_candidates(
	    network: BooleanNetwork,
	    candidates: list[dict[str, int]],
	    terminal_restriction_space: BinaryDecisionDiagram,
	    max_iterations: int,
	    ensure_subspace: dict[str, int] = {}
	) -> list[dict[str, int]]:

Currently, the `_preprocess_candidates` function is common for the two cases: 1) inside minimal trap spaces and 2) outside maximal trap spaces.
In `mtsNFVS`, these two cases are treated separately.

For Case 2, we expect to reach the best case where the candidate set is empty after finishing the `_preprocess_candidates` function.
The current code is OK.

	for state in candidates:
		state_bdd = state_to_bdd(state)
			...
	    for _ in range(max_iterations):

For Case 1, we should change a bit.

	for _ in range(max_iterations):
		for state in candidates:
			state_bdd = state_to_bdd(state)
			...
		
		if len(candidates) == 1:
			break

I will do this important task.

### The problem of secondary source nodes

Consider the Boolean network

	A, A
	B, B & A
	C, C & A

By running the working version, we obtained

	# attractors = 5
	# nodes = 11, time = 0.01

The number of SD nodes should be only 7.
The problem is here that when fixing node `A` to 1, node `B` and node `C` become source nodes, but `Trappist` now does not recognize them.

We need to deal with this issue.

### Function  `find_minimum_NFVS`  is not deterministic

I agree with Sam's suggestion <https://github.com/jcrozum/nfvs-motifs/issues/36>.


### Partial order reduction

> {Sam}: Listen, I've been playing with the succession diagram generator and I have a question: Have you considered incorporating something like partial order reduction to the structure of the diagram? Because right now, some models seem to contain a lot of combinatorial explosion that blows up the diagram quite a bit. Essentially, I am thinking about an approach similar to the way you treat input nodes in pystablemotifs, but more systematic (i.e. it could be applied at any level, as long as the stable motifs are sufficiently "independent"). My current intuition is that for attractor search it should work fine, just with smaller succession diagrams. But it may influence the control algorithm (I believe the algorithm could still work, albeit with some modifications).



It would be great if partial order reduction can work here.



### Finish some todo tasks

> TODO: properly create terminal restriction space; for now, just avoid stable motifs

> TODO (1): There are multiple places where we build the symbolic encoding, and it
>    will surely introduce extra overhead. We might want to just create the encoding
>    once and then pass it around.

> TODO (2): We could probably make this algorithm slighlty less random by doing
>    a limited version of symbolic reachability. I.e. instead of simulating just one
>    state transition in each step, compute the whole successor BDD and then test 
>    against that. Once the BDD becomes too large after several steps, we can just 
>    pick a single state from it and start again. Sam: I'll add a version of this
>    later, once we can actually benchmark how it performs :)
