## Algorithm

Here, we give a high-level description of the main attractor detection algorithm. 

*Sam: For now, I am using symbolic reachability in AEON to express some of the ideas. In the future, we can swap the reachability for a different technique, but if we want to have a quick prototype, we can just use it right away without almost any prerequisites.*

### Description

**Input:** A Boolean network $F$, a *negative-relaxed* network $F^{\oplus}$, a *trap space* $S \in \mathbb{B}_\star^n$.

> Sam: A "negative-relaxed" network is a variant of $F$ where all negative cycles in the signed influence graph of $F$ are eliminated. This can be done by fixing the values of variables that together form a (minimal) feedback vertex set of the signed influence graph of $F$ (the NFVS method by Van Giang). Alternatively, it can be also done on the level of fixing individual regulations within the update functions $f_i$ (we discussed this with Van Giang, but, afaik, it is not implemented). This can potentially lead to a less restrictive network with fewer fixed-points, since the variables do not necessarily have to be fixed to constant values, but can retain some (positive) influences. For now, the original NFVS variant should be sufficient though.

**Output:** A collection $[A_1, \ldots, A_k]$ of all *asynchronous attractors* within $S$.

**Alternative/extended outputs:** The "successors" of the space $S$ in the *succession diagram* of $F$. Also, instead of the full attractor sets $A_i \subseteq \mathbb{B}^n$, we can just output one *representative* state $A_i^{\odot} \in \mathbb{B}^n$ for each attractor, or a minimal space $A_i^{\square} \in \mathbb{B}_\star^n$ that contains the attractor (note that $A_i^{\square}$ is not necessarily a trap space—a minimal trap space $M$ can be larger than any $A_i^{\square}$ of attractors in $M$).

1. Value propagation (*percolation*): Compute $S'$ as the space where all fixed values in $S$ are propagated. This can be simply performed by abstract evaluation over Boolean expressions.

   > Sam: I have implemented this recently in AEON, should be easy to transfer into a new implementation.

   *Lemma/Conjecture:* There is no attractor in $S \setminus S'$, since every state in $S \setminus S'$ can reach a state in $S'$ in finitely many steps and $S'$ is a trap space.

   By the lemma above, we can recursively call the main procedure on $S'$ without considering any of the remaining states in $S$. 

2. Find a collection of maximal trap spaces within $S$. Let's denote them $M_1, \ldots, M_k$.

   *Lemma/Conjecture:* Note that maximal trap spaces can have non-empty intersections, in which case we might recursively cover the same sub-space multiple times. Alternatively, we can consider a sequence of spaces that are maximal without intersecting any other space in the sequence (this is slightly different than "global" maximality). Note that this definition does not guarantee that such sequence is unique, but that is not an issue: the method should work given any such valid sequence, as these "disjoint maximal" trap spaces still cover all minimal trap spaces within $S$ (this is in fact the claim of this lemma).

   In the following, I am thus silently assuming that $M_i$ are pair-wise disjoint. However, the algorithm also works if $M_i$ are simply the maximal sub-spaces within $S$.

   > Sam: At the moment, this should be performed by Trappist. As far as I understand, Trappist can now compute maximal trap spaces. It is my understanding that it could also be quite easily transformed to iteratively compute the "disjoint maximal" trap spaces if desired.

3. First, assume that $M_1, \ldots, M_k$ is empty. That is, there is no sub-space of $S$ that is also a trap space. Conversely, $S$ is a minimal trap space. 

   Now, compute "attractor candidate states" $N \subseteq S$ that correspond to the fixed-points within $F^{\oplus}$. These "cover" every real attractor in $F$. Then, use some abstract/symbolic reachability method to pre-order the elements of $N$ based on mutual reachability. The terminal elements of this pre-order are the "true" representatives of the individual attractors.

   > Sam: As far as I recall, NFVS by Van Giang uses Petri net unfolding to obtain this pre-order on $N$. To simplify the method, we may start by using symbolic BDD reachability first and only resort to other methods if this fails. In my experience, once the sub-space where the attractor resides is known, symbolic reachability is often quite fast (this is what we also use in AEON: we try to converge towards such subspace using reductions, but most of the time is spent on these reductions, not on detecting the attractor in the final subspace).
   >
   > Another alternative may be that if symbolic/unfolding analysis takes too long, we return the whole $S$ as a "reasonable approximation" of the resulting attractors.

   Finally, we return a collection of attractors $A_1, \ldots, A_k$ (or their representatives) that reside within this minimal trap space $S$.

4. Now, assume that $M_1, \ldots, M_k$ is not empty. In that case, recursively compute the attractors within each $M_i$, resulting in a sequence $A_1, \ldots, A_t$. These are the attractors that reside within the maximal trap spaces. However, we still have to compute attractors that are outside of the maximal trap spaces.

5. For this, again compute the set $N \subset S$ as the collection of fixed-points with respect to $F^{\oplus}$ such that $N \cap M_i = \emptyset$ for every $M_i$ (i.e. only compute candidate fixed points that are not within any of the maximal trap spaces).

   > Sam: My understanding is that any reasonable SAT/SMT/ASP solver should be able to compute $N$ reasonably fast, assuming the list of $M_i$ is not very long. However, if the list of $M_i$ is too long for the solver, it is probably too long for step (4) as well. So realistically, it shouldn't be an issue.

   Subsequently, we compute $N' \subseteq N$ consisting of states that cannot reach any $M_i$. For this $N'$, we perform the same procedure as in step (3). That is, we pre-order them based on reachability and identify the actual attractors (or attractor representatives).

   > Sam: Again, symbolic reachability could be used as the first "simple" method for computing $N'$. The advantage is that $S$ and all $M_i$ have a very concise symbolic representation and there is only a limited number of variables that can be updated by the reachability process (anything that is free in $S$ but fixed in some $M_i$). Hence the practical complexity isn't really the number of network variables, but rather the number of variables in which $S$ and $M_i$ differ. If we find anything faster or more reliable, we can of course use that instead.

   In this step, we obtain the "trap avoidant" attractors within $S$, let's denote them $A_{t+1}, \ldots, A_{k}$. Subsequently, we return the full sequence $A_1, \ldots, A_t, A_{t+1}, \ldots, A_k$. The algorithm is then complete.

##### Notes

In step (5), there are other possibilities of computing $N'$, for example based on the trap spaces of the time reversal network $F^{-}$. However, these are incomplete (*Sam: as far as I know*), that is, $N'$ cannot be computed fully solely based on $F^-$. We can thus later consider speeding up the computation of $N'$ using these time-reversed trap spaces, but initially this may not be necessary.

### Pseudocode

Note that the pseudocode does not track the scheme above exactly, but the general ideas should hopefully be clear.

```python
def attractors(F, F_plus, S):
  # (1) Check if values in S can be updated through
  # value propagation, and if so, restart search.
  P = propagate(F, S)
	if P != S:
		# There are no attractors in S \ P.
		return attractors(F, F_plus, P)
	
  # (4) Recursively process all maximal (or disjoint 
  # maximal) subspaces of S.
	trap_attractors = []
	traps = BDD.empty_set()
	for M in maximal_disjoint_subspaces(F, S):
		trap_attractors += attractors(F, F_plus, M)
		traps = BDD.union(traps, M)
	
  # (5a) Compute the "restriction space" in which no
  # trap avoidant attractor within S can reside.
  # For now, this is an exact method, but could be
  # approximate in later versions.
	traps_bwd = symbolic_backward_reachability(F, traps, S)
	
  # (5b) If the restriction space is the whole S,
  # we are done.
	if traps_bwd == S:
		# Every state in S can reach some M_i and thus
		# there are no trap avoidant attractors.
		return trap_attractors;
		
  # (5c) Compute the candidate set (N' in step (5))
  # based on `F_plus` and `traps_bwd`.
	N = restricted_fixed_points(F_plus, S, traps)
	N = N.minus(traps_bwd)
	
  # (5d/3a) Go through all the candidates and compute
  # the attractor for each candidate (if applicable),
  # plus also eliminate all candidates that can reach
  # this candidate (these cannot be in attractors).
	non_trap_attractors = []	
	while N.len() > 0:
    pivot = N.pick()
    bwd = symbolic_backward_reachability(F, pivot, S)
    # Set fwd actually does not have to be computed fully.
    # We can stop once it is clear pivot is not in an 
    # attractor (i.e. fwd is not a subset of bwd).
    fwd = symbolic_forward_reachability(F, pivot, S)
    # States in bwd cannot be attractor candidates.
    N = N.minus(bwd)
    if fwd.minus(bwd).is_empty():
      non_trap_attractors.append(fwd)
      
	return trap_attractors + non_trap_attractors
  	
# Value propagation of `S` with respect to the update
# rules of the given `network`.
def propagate(network , S):
  # Should be easy to implement manually using abstract
  # evaluation over Boolean expressions.

# Sequence of maximal disjoint subspaces in `S`.
def maximal_disjoint_subspaces(network, S):
	# Implemented based on Trappist.

# Symbolic set of all states that can reach `pivot` within 
# the given `space`.
def symbolic_backward_reachability(network, pivot, space):
	# Implemented based on AEON.

# Symbolic set of all states reachable from `pivot`
# within the given `space`.
def symbolic_forward_reachability(network, pivot, space):
	# Implemented based on AEON. 

# All fixed-points of the given `network` that are
# within the given `pos_space`, but outside of 
# the given list of `neg_spaces`.
def restricted_fixed_points(
  network, pos_space, neg_spaces
):
	# Implemented using any reasonable SAT/SMT/ASP solver.
```





