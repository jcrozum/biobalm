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

We can also consider the *time reversal* of $F$, denoted $F^-$. The update functions $f_i^-$ are defined such that $f_i^-(x) = \neg f_i(x_1, \ldots,x_{i-1},\neg x_i,x_{i+1},\ldots x_n)$. The result is a network with a reversed state-transition graph compared to $F$. That is, $x \to y$ in $STG(F)$ if and only if $y \to x$ in $STG(F^-)$.

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

##### Single-node driver set $\Delta$

For a set of maximal trap spaces, $\mathscr{M}$, $\Delta=\Delta(\mathscr{M})$ is the set of all variable-value pairs $(x_i,s_i)$ such that there exists a maximal trap space $M\in\mathscr{M}$  that contains the subspace obtained by percolating $x_i=s_i$.

Rozum et al. 2021 proved that if an attractor exists outside of the set $\bigcup_{M\in\mathscr{M}} M$, then that attractor must lie within the subspace defined by $\neg \Delta$, i.e., each state in the attractor must satisfy $x_i=\neg s_i \space \forall(x_i,s_i)\in\Delta$.

> Giang: I have two questions with this definition.

First, let us consider the BN shown in Figure 4 of Rozum et al. 2021.
$$
f_A = (\neg A \land \neg B) \lor C\\

f_B = (\neg A \land \neg B) \lor C\\

f_C = A \land B
$$
This BN has one max. trap space $\{A = 1, B = 1, C = 1\}$.

In the paper, $\Delta = \{C = 1\}$ (i.e., $\neg \Delta = \{C = 0\}$) and we have $R(X) = \neg C \land (\neg A \lor \neg B)$.

$R(X) = 1$ contains three states 000, 010, and 100.

However, from the definition, $\Delta$ can be $\{A = 1, C = 1\}$ (i.e., $\neg \Delta = \{A = 0, C = 0\}$).
> Jordan: (A,1) is not an element of $\Delta$ according to the definition because percolating $A=1$ gives the subspace $1\star\star$, which does not lie within the subspace $111$. However, percolating $C=1$ gives $111$, so $(C,1)\in\Delta$.

We have $LDOI(A = 0) = \{C = 0\}$.

Then, we get $R(X) = (\neg C \land (\neg A \lor \neg B)) \land (\neg C \land \neg A \land B)$.

$R(X) = 1$ contains only one state 010.
> Jordan: With the correct $\Delta$, $R(X)=1$ on the set $\{(A,B,C): C=0, A\wedge B=0\}$, i.e., the set $\{000,010,100\}$, which is the full motif-avoidant attractor in this case (but may contain additional states in general).

The BN under the fully asynchronous update has one motif-avoidant attractor including three states 000, 010, and 100.

Hence, the statement that "All motif-avoidant attractor states satisfy $R(X) = 1$" seems not hold.

Jordan, could you please check this?



Second, it seems that $\Delta$ should be as large as possible to reduce the number of states covered by $R(X) = 1$ (also the terminal restriction space).

How does pystablemotifs compute $\Delta$?



> Giang: Sorry for my misunderstanding. But I think the definition may be that "For a set of maximal trap spaces, $\mathscr{M}$, $\Delta=\Delta(\mathscr{M})$ is the set of all variable-value pairs $(x_i,s_i)$ such that **for every $(x_i,s_i)$** there exists a maximal trap space $M\in\mathscr{M}$  that contains the subspace obtained by percolating $x_i=s_i$." With this definition, $\Delta$ contains **all** such $(x_i,s_i)$, i.e., it is as large as possible.
> > Jordan: I like this change, and the definitions are equivalent. I think the confusion came from whether the "such that" was inside or outside the set comprehension. Your definition eliminates this confusion.
> >
> > > Giang: Now, I completely understand $\Delta$. Thank you.
>
>
> Giang: In addition, I think the computation of $\Delta$ is not hard. I guess that you consider all `n` possible single-node driver sets, for each you percolate and check if it is contained in a max. trap space.
> > Jordan:  Yes, this is basically correct, but it's actually $2n$ possible single-node driver sets (because there are two states per variable).
> >
> > > Giang: My bad. It should be $2n$ indeed.



> Giang: Once we have obtained $\Delta$, we can get $R(X) = 1$ in which all motif-avoidant attractors must reside. I guess that you are trying to find a **finer** formulation for $R(X)$ based on time-reversal. This means to prune more parts of the state space. Is it right?

##### Minimum driver node of a stable motif

*Driver Set*: Given a maximal trap space $M$, a *driver set* of $M$ is any set $S$ of variable-value pairs $(x_i,s_i)$ such that the subspace obtained by percolating $S$ lies within the subspace $M$.

> Giang: I think there is another possible definition that is "Given a maximal trap space $M$, a *driver set* of $M$ is any set $S$ of variable-value pairs $(x_i, s_i)$ such that **all minimal trap spaces** of $N_S$ lies within the subspace $M$ where $N_S$ is the BN obtained from the original BN by fixing all node $x_i$ to $s_i$." This definition is stronger than the old one. Indeed, if a driver set satisfies the old definition, it also satisfies the new definition, whereas if a driver set does not satisfy the old definition, it may satisfy the new definition. For example, consider the example BN shown in Slide "Improve the accuracy" of my presentation. Following the old definition, we should get the minimum control policy $\{x_1 = 1, x_3 = 0\}$. Following the new definition, we should get the minimum control policy $\{x_3 = 0\}$ that does not satisfy the old definition. How do you think about the new definition?

*Internal Driver Set*: A driver set $S$ of a maximal trap space $M$ is an *interal driver set* if and only if the subspace defined by $S$ contains the subspace $M$. Equivalently, a driver set $S$ is an internal driver set if every variable-value pair in $S$ corresponds to a fixed variable of $M$.

*Minimal Driver Set*: A driver set $S$ of a maximal trap space $M$ is minimal if there does not exist a subset $T$ of $S$ that is also a driver set of $M$.

There are thus a few versions of the problem of finding driver nodes for a maximal trap space:
1. Find all minimal driver sets
2. Find all minimal internal driver sets
3. Find all smallest driver sets (i.e., the minimal driver sets with the fewest number of variable-value pairs)
4. Find the smallest internal driver sets
5. Find any minimal driver set

By default, `pystablemotifs` tries to solve problem 2. There is an option to have it attempt problems 3 or 5.



> Giang: Regarding the target control process, I assume that we want to find all minimum control policies that drive any initial state of the network to a desirable min. trap space $m$. I think the benefit when considering the sequences of max. trap spaces leading to $m$ and only the internal nodes (i.e., the history internal control) is that the number of possible solutions is only $O(\sum_{i = 1}^k (2^{|D(M_i)|}))$ instead of $O(2^n)$ where $D(M_i)$ is the set of fixed variables of trap space $M_i$, $n$ is the number of nodes of the original BN, and $\sum_{i = 1}^k |D(M_i)| \leq n$ (leading to $\sum_{i = 1}^k (2^{|D(M_i)|}) << 2^n$). However, if we consider both internal and external nodes with respect to a max. trap space, this benefit disappears as the number of possible solutions is $O(2^n)$. Hence, I think if we consider both internal and external nodes, we should consider **directly** the target min. trap space $m$ instead of the sequences of max. trap spaces leading to it. This would be a better approach.

> Giang: Jordan, could you please add the formal definition of this problem?



##### TBD

Regulatory/signed influence graphs, stable motifs, feedback vertex sets,...

[bryant-86] Bryant, Randal E.. “Graph-Based Algorithms for Boolean Function Manipulation.” *IEEE Transactions on Computers* C-35 (1986): 677-691.