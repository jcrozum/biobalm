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
