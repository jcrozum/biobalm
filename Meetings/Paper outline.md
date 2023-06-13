# Trap spaces vs. attractors: interplay of positive and negative feedback in asynchronous Boolean dynamics

## Motivation

Boolean networks are one of the fundamental modelling frameworks in systems biology. However, the theoretical as well as practical aspects of update concurrency within Boolean networks is still a widely discussed topic. Here, an important distinction is between attractors and minimal trap spaces. 

An **attractor** is a terminal (bottom) SCC of a state-transition graph induced by a particular update mode. Meanwhile, a **trap space** is an inclusion-minimal hyper-cube that cannot be escaped under *any* update mode. As such, in any update mode, a trap space is *guaranteed* to contain at least one attractor. However, (1) multiple attractors can appear within a single trap space, and (2) attractors can sometimes appear outside of a minimal trap space.

It is known that in practice, minimal trap spaces approximate attractors in most update modes *very well* (and they are known to be *exactly* attractors under the *most permissive* update). However, it is unclear exactly how "good" is this approximation, how "rare" are *trap-avoidant* attractors, and how this property "scales" across different asynchronous update modes. Furthermore, it is unclear whether "biologically motivated" restrictions of the random BN framework such as local-monotonicity of update functions impact the presence of such attractors.

## Scope

We consider four types of update concurrency (in the order of "permissiveness"):

- "Plain" asynchronous update.

  > "Plain" asynchronous networks can contain motif-avoidant attractors:
  >
  > ```
  > f_a = a <=> b
  > f_b = a <=> b
  > ```

- "Generalised" asynchronous (also "set" asynchronous): Any combination of variables can update together.

  > "Generalised" asynchronous networks can contain motif-avoidant attractors:
  >
  > ```
  > f_a = (!a & !b & !c) | (a & b & c)
  > f_b = (!a & !b) | (a & b & c)
  > f_c = (a & !b & !c) | (c & a) | (c & b)
  > ```
  >
  > Should have a 111 fixed point and the remaining states are one big SCC. This actually also works as "plain" asynchronous, but is ugly... Does anyone know about a 2D network under gen. async. updates? I don't think it's possible, but a 2D example would be nicer.

- "Interval" asynchronous update: New value becomes "visible" with a non-deterministic delay, but is locked in place once the update happens (i.e. an update cannot stay "invisible"):

  - $00 \to 10$ / $11 \to 01$ based on update function: the function changes the "stored" value, but only if stored = visible.
  - $10 \to 11$ / $01 \to 00$ stored value propagates to visible value any time.

  > Boolean Networks: Beyond Generalized Asynchronicity; Thomas Chatain, Stefan Haar, and Loic Pauleve
  >
  > Authors shows that "interval" admits every transition available in the generalised asynchronous case. But does not show anything about the inverse direction.
  >
  > The example above does not work, because I can do:
  > $$
  > 0,0,0 \to^{f_b} 0,0(1),0 \to^{f_a} 0(1),0(1),0 \to^{\epsilon} 1,0(1),0 \to^{f_c} 1,0(1),0(1) \to^{\epsilon} 1,0(1),1 \to^{\epsilon} 1,1,1
  > $$
  > This shows that the "interval" dynamics contains transitions that are not present in the generalised asynchronous one. Furthermore, it clearly also admits fewer edges than "linear cuts" because in the "linear cuts" case, I can do $0,0(1),0 \to 0,0,0$, but I cannot do this in the "interval" case. The question is whether this is useful for something. 

- Asynchronous update with "linear cuts": New value becomes "visible" with a non-deterministic delay, but it can also "fall back" to the previous value:

  - $00 \to 10$ / $10 \to 00$ / $01 \to 11$ / $11 \to 01$ based on update function: the function can change the "stored" value any time, regardless of what is observed.
  - $10 \to 11$ / $01 \to 00$ stored value propagates to visible value any time.

  > Linear cuts in Boolean networks Aurelien Naldi , Adrien Richard, and Elisa Tonello
  >
  > In the paper, the concept is slightly different: They assume the variables facilitating the "linear cut" are already present in the network. However, one can easily imagine just adding a "cut variable" for every network variable and the result would be the dynamics described above.
  >
  > At the same time, the authors show that "linearly cuttable" network has no motif-avoidant attractors and that every minimal trap space contains exactly one attractor. However, this attractor can be smaller than the trap space (as opposed to the "most permissive" update).
  >
  > Furthermore, linear cuts over-approximate "single threshold refinements" but not general multi-valued dynamics.

- "Most permissive" update: A value change goes through an intermediate state $\star$ that essentially means "any value". The important distinction is that such $\star$ can be perceived as 0/1 in different contexts (previous modes had exactly one visible value at any time):

  - $0 \to \star$ / $1 \to \star$ based on update function ($\star$ can be seen both as 0 and 1 in other update functions).
  - $\star \to 0$ / $\star \to 1$ any time.

  > The reason why this admits even more behaviour is that the $\star$ symbol can mean anything to anyone, so it can even represent both 0 and 1 in the same function. 

> In theory, the "linear cut" and "interval" update modes can be encoded into a "plain" asynchronous network of $2n$ variables (without a blow up in the influence graph) and analysed that way. For "most permissive" update, we probably don't need to do anything special because there we just care about trap spaces (but it can be also encoded into $2n$ variables, it just requires more involved manipulation of update functions).

Aside from the "hierarchy" above, I am not aware of any work linking locally-monotonic functions to motif-avoidant attractors. It **seems** that no one has an example of a motif avoidant attractor based on a locally-monotonic network, but it just might be because high dimensionality is needed to produce such attractor?

> I just realised I can probably express this as a HCTL property and run this as an exhaustive model checking query on fully unknown network with a complete influence graph. This probably won't scale beyond a few variables, but if there is a small counterexample, this should find it. I'll report back with results :)

We study the following network types:

- "Real world" networks from the BBM benchmark. For networks with inputs, we either test all input valuations, or sample a reasonably high number of valuations randomly (when there are too many inputs).
- Four "classes" of random N-K networks with $K=2$:
  - Exactly $K=2$ inputs on each variable vs. Poisson distribution of degrees, but $K=2$ average **(do we allow inputs?)**.
  - Random Boolean functions vs. locally-monotonic functions **(essential inputs?)**

Ultimately, the questions are:

- Can interval networks contain motif-avoidant attractors? If yes, what is the "scaling" of such attractors compared to "plain" asynchronous?
- Can locally-monotonic networks contain motif-avoidant attractors? If yes, what is the "scaling" of such attractors compared to arbitrary networks?
- Is the attractor scaling the same for $K=2$ but different in-degree distribution?

> I know that this is not necessarily a contribution, but it would be nice to have a single paper where these effects on attractor scaling and trap-avoidant attractors are summarised coherently instead of having it scattered across five papers.

## Contributions

1. We formalise the concept of **partial succession diagrams** (name pending :)).

2. We introduce an attractor detection method based on partial succession diagrams, negative FVSes and symbolic methods. We show that this method performs better than other available methods.

   1. The method starts with a "minimal viable" succession diagram that only guarantees a single path to each minimal trap space.
   2. Then, based on the candidate sets given by the FVS method, it can further expand the succession diagram to *separate* as many candidates into smaller sub-spaces as possible.
   3. Then, once trap spaces can no long separate the candidates, we use heuristics/static analysis (simulation, pint) and exact methods (symbolic reachability, TGR, mole) to prune the candidates.
   4. Optionally, symbolic reachability can be used to compute the exact attractor set.

   > This is just a template trying to separate the algorithm into distinct steps for easier presentation. The final code can be different. Also, other optimisations can be part of the algorithm depending on what we measure performs best, this is just a broad outline. The core is the idea that we start from a "minimal viable" diagram, expand it until it no longer "helps" with attractor checking and then run other methods on the resulting decomposed state space.

3. For real-world networks, we study the following questions:

   1. **What are the differences between min. trap spaces and attractors within them?** Is there more than one attractor? Is there a smaller sub-space that is not a trap but contains the attractor? How many states of the trap space are attractor states?
   2. **Are there any trap-avoidant attractors?** If so, how many? What is the distance of such attractors to minimal trap spaces? (in terms of succession-diagram "steps" or hamming distance to the minimal trap spaces that contains the attractor).

   Depending on what we can actually prove, we might want to test this for both "plain" and "interval" update.

4. We do the same thing for different flavours of random networks, but here the goal is to explore the "scaling" of such properties with increasing network size.

5. We should explore the scaling of average succession diagram size with network size? The absolute node count probably grows very quickly, but the **height** of the succession diagram should be very small even for large networks (it is clearly always $\leq n$). It could be interesting to figure out if it is closer to $log(n)$ or $\sqrt{n}$ or something entirely different.

#### Other contributions

1. It seems slightly unlikely but possible that some of these scenarios will not allow trap-avoidant attractors at all (e.g. no trap avoidant attractors with interval update). It would be very cool if we could **prove** something like that formally. 
2. **Control:** We could design a succession diagram expansion scheme that is specifically targeting control. My understanding is that we basically want to expand the succession diagram to the largest trap spaces that are still within the *strong basin* of our target attractor/phenotype in order to minimise the size of perturbation. This should be achievable. However, I am not sure how to reconcile this with the rest of the proposed contributions. Maybe something like **"Does presence of motif-avoidant attractors hinder control?"**
3. General properties of attractor/trap space scaling in different random models. Different papers explore different random classes of models, but AFAIK there isn't a paper that would compare these somewhat sensibly. (see one of the notes above)

## TODO list (work in progress)

- Paper (some of this may end up as supplementary material):
  - [ ] Definitions of trap space (minimal/maximal) and attractor (of a generic "STG").
  - [ ] Define value percolation (w.r.t. to spaces).
  - [ ] Define succession diagram and partial succession diagram (some nodes are unexpanded "stubs").
  - [ ] Algorithm for building a "minimal viable" partial succession diagram using trappist (i.e. one path to each minimal trap space and nothing more).
  - [ ] Algorithm for "NFVS guided expansion" of the succession diagram for attractor detection.
  - [ ] Describe how attractors are detected once the succession diagram cannot be expanded further (simulation, reachability, transition guided reduction, etc.) and other algorithmic stuff. *This is still rather open because we are experimenting with best/fastest way of doing this, so this is mostly a ""task template" for later.*
  - [ ] Define the different "flavours" of asynchronous update that we are interested in. Maybe also define the translations to a Boolean network of size $2n$ (somebody will ask for it in the supplement for sure).
  - [ ] Define the classes of random models we are interested in.
  - [ ] Try to find examples of "interval" and "linear cut" networks with motif-avoidant attractors, ideally with locally-monotonic update functions. (If we can't find any, this just got more interesting)
  - [ ] To be continued based on experiments...
- Implementation:
  - [ ] Succession diagram expansion based on candidates obtained through NFVS method. Test how large are the succession diagrams with this? (Will be probably larger than with symbolic pruning, the question is how much larger?)
  - [ ] Build/obtain a generator for random networks that we can tune sufficiently for our needs.
  - [ ] Implement the translation for the two asynchronous dynamics to "plain" asynchronous.
  - [ ] We will for sure need to visualise succession diagrams at some point. It would be super useful if we could export a succession diagram to a `.dot` file or directly into a Tikz figure.
  - [ ] Start experimenting with control...
  - [ ] [Sam] Try replacing symbolic reachability with transition-guided reduction.
  - [ ] [Sam] Add symbolic percolation directly into AEON.py.
  - [ ] [Sam] A nicer Petri net translation based on AEON BDDs and a purpose-made data structure.
- Experiments and benchmarks:
  - [ ] Build the benchmark sets of random networks (once the generator is available).
  - [ ] Build the benchmark set of real-world networks (this mainly includes sampling different input valuations).
  - [ ] Run attractor detection on the benchmark sets (if we run into performance problems, we should re-explore the main algorithm). This should yield upper bound on SD size, number of motif-avoidant attractors, number of attractors in each trap space, and ideally also the symbolic set for each attractor (at least for smaller networks). *This step is only concerned with the results, not performance.*
  - [ ] Benchmark attractor detection on standardised hardware using the final implementation, comparing to AEON, NFVS and pystablemotifs. 