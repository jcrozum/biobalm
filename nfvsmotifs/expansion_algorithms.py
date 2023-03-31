from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from biodivine_aeon import BooleanNetwork, SymbolicAsyncGraph, VariableId, ColoredVertexSet # type: ignore

from nfvsmotifs.SuccessionDiagram import SuccessionDiagram
from nfvsmotifs.space_utils import symbolic_space
from nfvsmotifs.trappist_core import trappist

DEBUG = False
SYMBOLIC_PRUNING = True

def trap(
    network: BooleanNetwork,
    stg: SymbolicAsyncGraph,
    variables: list[VariableId],
    initial: ColoredVertexSet,
) -> ColoredVertexSet:
    """
    This method takes an `initial` symbolic set of states and iteratively eliminates all states from 
    which the set can be escaped using transitions utilizing one of the given `variables`.

    The method can terminate prematurely if the itermediate set becomes too large (in terms of BDD 
    nodes). In such case, the result is a subset of the `initial` set, but a superset of the 
    "optimal" result. Currently, the cut off is set to `100_000 * network.num_vars()` (i.e. 
    roughly ~100kB of RAM used per one network variable). Right now, there is no limit on the
    number of iterations, but as long as the BDD is sufficiently small, it's usually beneficial
    to just let the method run as long as necessary.

    For optimal performance, the `variables` list should be ordered in accordance to the internal
    variable ordering of the associated BooleanNetwork, but this is not strictly required.
    """
    # Currently, `set.symbolic_size()` should report the approx. size of the BDD in bytes,
    # so this should literally be 100kB/variable.
    size_limit = 100_000 * network.num_vars()
    result = initial

    keep_working = True
    while keep_working:
        keep_working = False
    
        for var in reversed(variables):
            can_post_out = stg.var_can_post_out(var, result)

            if not can_post_out.is_empty():
                keep_working = True
                result = result.minus(can_post_out)
                
                if DEBUG and result.symbolic_size() > 100_000:
                    # If the BDD is overall larger than 100kB we should probably
                    # start logging progress.
                    print(f"BDD node count: {int(result.symbolic_size() / 10)}; Remaining state space: {result.cardinality()}.")                

                if result.symbolic_size() > size_limit:
                    if DEBUG:
                        print("Symbolic reachability terminated early.")
                    return result
                
                break

    return result

def simplified_dfs_expansion(sd: SuccessionDiagram):
    """
    This algorithm expands the entire succession diagram in a DFS manner, but it tries to skip as many nodes
    as possible to reduce the size of the resulting diagram.

    In particular, the resulting diagram satisfies that:
     - Every minimal trap space is still present.
     - Every node is either a "stub" node with no successors, or is expanded in the same way as in the "full" diagram.
     - Following inductively from the previous property, every path leading to a minimal trap space is also present 
     in the "full" diagram.

    Furthermore, as a byproduct of the expansion procedure, we can often prove that there are no motif-avoidant
    attractors in the sub-space of the currently expanded node. As such, if you use this expansion procedure,
    it can greatly simplify any subsequent attractor search.

    The main idea of the procedure is the following: For each expanded node, we only need to expand stable motifs
    until all states corresponding to the node are either (a) direct members of one of the expanded stable 
    motifs, or (b) can reach one of the expanded stable motifs. Any additional stable motifs are redundant, as
    they would only re-explore the same set of states as already covered, or one of the states where we've proved
    it isn't possible to discover a new minimal trap space or attractor. Further details regarding correctness
    are present within the algorithm.
    """
    root = sd.root()
    root_space = sd.node_space(root)   # This is after percolation, so it might not be the whole state space.
    root_space_symbolic = symbolic_space(sd.symbolic, root_space)

    if len(root_space) == sd.network.num_vars():
        # A special case where the network has one fixed-point that is
        # resolvable through percolation. Not super common but happens sometimes
        # and will cause the main algorithm to do weird things...
        sd.expanded.add(root)
        return

    root_motifs = trappist(sd.petri_net, problem="max", ensure_subspace=root_space)
    # In all instaces, we want to have the "smallest" stable motif (in terms of state space; 
    # i.e. the most fixed variables) as the last entry of the list, because expanding it should 
    # lead us to the minimal trap spaces as fast as possible. But this not always optimal.
    # In any case, we should sort the results *somehow* in case the results from clingo are 
    # sorted non-deterministically.
    root_motifs = sorted(root_motifs, key=lambda x: len(x))

    # Each stack entry holds three items:
    # - ID of the node that is being expanded.
    # - List of stable motifs that still need to be expanded.
    # - Symbolic subset of the node space that hasn't been covered so far by any (expanded) child.
    # Here, "covered" means that either (a) it is a state inside the child's stable motif, or (b) it is a state
    # that we have proved can reach the child's stable motif.
    stack: list[tuple[int, list[dict[str, int]], ColoredVertexSet]] = [(root, root_motifs, root_space_symbolic)]

    # Note that in a normal DFS, we would also have to check if the node is already visited, but not 
    # fully expanded (i.e. "on stack"). However, we don't have to do it here, because we know that the 
    # succession diagram is always acyclic. So every node we "re-discover" is either completely fresh,
    # or already fully expanded.

    while len(stack) > 0:
        node, motifs, remaining = stack[-1]
        node_space = sd.node_space(node)

        if DEBUG:
            print(f"[{len(stack)}] Process node {node} with {len(motifs)} motifs and {remaining.cardinality()} states.")

        if len(motifs) == 0:
            sd.expanded.add(node)
            stack.pop()

            if remaining.is_empty():
                # If the set of remaining states is empty, this node cannot contain any motif-avoidant attractors.
                # Hence we can mark it as "finished" for attractor search as well.
                sd.attr_expanded.add(node)

            if DEBUG:
                print(" >> Done. No more motifs to process.")
                if not remaining.is_empty():                    
                    print(f" >> Found trap avoidant candidate states: {remaining.cardinality()}")
            
            continue
        
        if remaining.is_empty():
            # This node has no remaining uncovered states, but still has unexplored stable motifs.
            # These stable motifs can safely become stubs.
            for motif in motifs:
                child = sd.ensure_node(node, motif)
                if child not in sd.expanded:
                    # It is possible that this process will re-discover an existing
                    # expanded node. In such case, it would be invalid to mark it as a stub.
                    sd.node_is_stub(child, True)
            if DEBUG:
                print(f" >> Done. {len(motifs)} stubs created.")

            # As before, there can be no attractors in this space because the remaining set is empty.
            sd.attr_expanded.add(node)
            sd.expanded.add(node)
            stack.pop()            
            continue

        # Pick the next motif to process and expand it.
        motif = motifs.pop()
        motif_symbolic = symbolic_space(sd.symbolic, motif)
        child = sd.ensure_node(node, motif)        
        child_space = sd.node_space(child)  # Again, this is different from `motif` because motif is un-percolated.
        child_space_symbolic = symbolic_space(sd.symbolic, child_space)            

        # If this new node only *percolates* to states that are already covered,
        # we don't have to expand it further.        
        # However, we should still remove the *unpercolated* motif from the remaining
        # set of states, as this set *can* actually interset with the remaining states.
        is_redundant = child_space_symbolic.intersect(remaining).is_empty()
        
        remaining = remaining.minus(motif_symbolic)
        if SYMBOLIC_PRUNING:
            # Only variables that are *not fixed* within the parent space are relevant.
            free_variables = [ x for x in sd.network.variables() if sd.network.get_variable_name(x) not in node_space ]
            remaining = trap(sd.network, sd.symbolic, free_variables, remaining)
        
        # Now that `remaining` is updated, we can update the value on stack.
        stack[-1] = (node, motifs, remaining)        

        if DEBUG:
            print(f" >> Expanding child motif into node {child}. Reduced remaining states: {remaining.cardinality()}.")

        if is_redundant:
            if child not in sd.expanded:            
                sd.node_is_stub(child, True)
                if DEBUG:
                    print(f" >> Done. Discovered a new child stub.")
            else:
                if DEBUG:
                    print(f" >> Done. Child is redundant, but already expanded.")
            continue
        
        if child not in sd.expanded:
            if len(child_space) == sd.network.num_vars():
                if DEBUG:
                    print(" >> Done. Child is a fixed-point.")
                sd.expanded.add(child)
            else:
                # Finally, if the child node is not redundant nor expanded, we can add it as a new
                # stack entry for expansion, using its own stable motifs.
                child_motifs = trappist(sd.petri_net, problem="max", ensure_subspace=child_space)
                if len(child_motifs) == 0:
                    if DEBUG:
                        print(f" >> Done. Child is a minimal trap space ({len(child_space)}/{sd.network.num_vars()} variables fixed).")
                    sd.expanded.add(child)
                else:                
                    child_motifs = sorted(child_motifs, key=lambda x: len(x))
                    # In the child node, we only need to cover the *percolated* state space, because
                    # we know there are no attractors/min. traps in the *unpercolated* state space.
                    stack.append((child, child_motifs, child_space_symbolic))
        else:
            # In this case, the only point of the child is to eliminate states in the `remaining` set
            # and we don't have to explore it further, because it has already been explored.
            if DEBUG:
                print(f" >> Done. Child is not redundant, but already expanded.")
