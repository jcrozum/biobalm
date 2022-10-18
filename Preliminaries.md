## Preliminaries

Just a place where to put basic definitions and assumptions so that we can reference them quickly when necessary...

##### Boolean networks

Let $\mathbb{B} = \{0, 1\}$. In the following, we are assuming a Boolean network $F$ over $n$ variables, consisting of $n$ Boolean update functions $F = \{ f_1, \ldots, f_n \}$, where $f_i: \mathbb{B}^n \to \mathbb{B}$. We then call $\mathbb{B}^n$ the *state space* of the network $F$, and any vector $x \in \mathbb{B}^n$ a state.

##### Boolean network representation

Note that there are different ways the network $F$ can be actually represented. The most common case is by using Boolean expressions (typically using $\neg$, $\land$, and $\lor$, but other operators like $\Rightarrow$, $\Leftrightarrow$, and $\oplus$ are possible as well). Other option is to use binary decision diagrams *[bryant-86]* that encode the function into a directed acyclic graph. 

The advantage of using expressions is the simplicity of the representation. The disadvantage is that the representation is not canonical, hence some operations may have non-trivial complexity. Some tools mitigate this to some extent by relying on normal forms like DNF/CNF, or lists of prime implicants, but these can sometimes be costly to maintain.

The advantage of BDDs is that many operations are easy, and the BDD is always canonical, assuming the BDD can be actually constructed. However, some very complex networks may contain functions that cannot be trivially translated to BDDs, while still being written using Boolean expressions.

>  Conclusion: For now we concentrate on Boolean expressions, and if there is a specific need for other means of representation, we can try, e.g., BDDs.

##### Asynchronous update

The evolution of the state of network variables is determined by the individual update functions $f_i$. We then consider that each function can be applied *asynchronously* (although, many results regarding trap spaces hold regardless of the update scheme). As such, we consider a state-transition graph $STG(F) = (\mathbb{B}^n, \to)$, where the vertices are the network states, and the transitions update individual variables using the functions from $F$. Formally, $x \to y \Leftrightarrow (x \not= y \land \exists i. y = x[i \gets f_i(x)])$ (here, $x[i \gets b]$ denotes the substitution of the $i$-th component of $x$ for the value $b$).

Note that the number of transitions in $STG(F)$ can be often trivially derived from the properties of $f_i$. In particular, if $f_i$ does not depend on the value of $x_i$, then exactly 1/2 of all the network states can perform a transition changing the value of $x_i$. 

The bottom (terminal) strongly connected components of $STG(F)$ are called *attractors* of $F$ and are the main focus of our investigation. Additionally, we use the term *weak basin* of attractor $A$ to refer to all states that can reach attractor $A$, and the term *strong basin* to denote the states that can only reach attractor $A$ (and no other attractor).

##### Time reversal

We can also consider the *time reversal* of $F$, denoted $F^-$. The update functions $f_i^-$ are defined such that $f_i^-(x) = \neg f_i(\neg x_1, \ldots, \neg x_n)$. The result is a network with a reversed state-transition graph compared to $F$. That is, $x \to y$ in $STG(F)$ if and only if $y \to x$ in $STG(F^-)$.

The attractors of the time reversed network $F^-$ are the *repellers* of $F$ (also, top SCCs), and other graph-related terms are reversed similarly.

##### Trap spaces

Let $\mathbb{B}_\star = \{0, 1, \star \}$. A *space*, or a *hypercube* in $\mathbb{B}^n$ is a vector $S \in \mathbb{B}^n_\star$, such that a state $x \in \mathbb{B}^n$ belongs into this space (written $x \in S$) if and only if $x_i = S_i$ for all $i$ where $S_i \in \{0,1\}$. In other words, $S_i = \star$ indicates that the variable is *free*, while $S_i \in \{0,1\}$ indicates that the variable has a prescribed *fixed* value.

A space $S \in \mathbb{B}^n_\star$ is then a *trap space* if for every state $x \in S$, and every $y \in \mathbb{B}^n$ such that $x \to y$ in the $STG(F)$, we have that $y \in S$. In other words, no state in $S$ can escape from $S$ using the transitions in $STG(F)$ (however, note that this property holds universally for all common BN updating schemes).

Observe that spaces form a partial order based on set inclusion. We can thus define a *minimal* trap space $M$ as the trap space that has no sub-space $M' \subset M$ where $M'$ is also a trap space. Similarly, we can define a *maximal* trap space as a space $M$ that has no super-space $M \subset M'$ such that $M'$ is also a trap space, and $M' \not= \star^n$. The last condition is necessary to obtain meaningful results, as $\star^n$ is trivially always a trap space that is a super-space of every other space.

Alternatively, we can define a minimal and maximal trap space *within* some predefined space $S$. In such a case, we consider $M$ to be minimal (or maximal) among the sub-spaces of $S$, but not equal to $S$.

##### Trap spaces and $STG(F)$

Note that (globally) minimal trap spaces have a strong correspondence with network attractors: every trap space contains at least one attractor. However, not every attractor must be contained in a minimal trap space. There can be asynchronous attractors that reside in larger trap spaces of the network.

Observe that a trap space $S$ in the time reversed network $F^-$ corresponds to a set of states for which there are no transitions entering $S$ (in $STG(F)$). That is, once $S$ is escaped, it cannot be entered again. However, note that this does not mean that attractors of $F$ cannot be contained in any trap space of $F^-$. As long as the attractor and the weak basin of the attractor are within the same (time reversed) trap space, there is no contradiction.

Finally, note that trap spaces of both $F$ and $F^-$ are SCC-closed with respect to $STG(F)$. That is, given a trap space $S$, every SCC $X$ of $STG(F)$ is either fully in $S$ (i.e. $X \cap S = S$), or outside of $S$ (i.e. $X \cap S = \emptyset$).

##### TBD

Regulatory/signed influence graphs, stable motifs, feedback vertex sets,...

[bryant-86] Bryant, Randal E.. “Graph-Based Algorithms for Boolean Function Manipulation.” *IEEE Transactions on Computers* C-35 (1986): 677-691.