from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from biodivine_aeon import SymbolicAsyncGraph, VariableId, ColoredVertexSet # type: ignore

from nfvsmotifs.SuccessionDiagram import SuccessionDiagram
from nfvsmotifs.space_utils import percolate_space, symbolic_space
from nfvsmotifs.trappist_core import trappist

DEBUG = False
SYMBOLIC_PRUNING = True

def trap(
    stg: SymbolicAsyncGraph,
    variables: list[VariableId],
    initial: ColoredVertexSet,
) -> ColoredVertexSet:
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
                    print(f"BDD node count: {int(result.symbolic_size() / 10)}; Reached {result.cardinality()}")                

                if result.symbolic_size() > 1_000_000:
                    return result
                break
    return result

def reach_bwd(
    stg: SymbolicAsyncGraph,
    variables: list[VariableId],
    initial: ColoredVertexSet,
    universe: ColoredVertexSet,
) -> ColoredVertexSet:
    result = initial    

    keep_working = True
    while keep_working:
        keep_working = False

        for var in reversed(variables):
            pre = stg.var_pre(var, result)

            if not pre.is_subset(result):
                result = result.union(pre)
                if DEBUG and result.symbolic_size() > 1_000_000:                    
                    print(f"BDD node count: {int(result.symbolic_size() / 10)}; Reached {result.cardinality()}/{universe.cardinality()}")
                keep_working = True

                # TODO: This limit should by somehow proportional to the size of the
                # network... something like 100_000 * variable_count seems reasonable?
                # Furthermore, we might want to have a limit on the total number of iterations.
                if result.symbolic_size() > 1_000_000:
                    return result
                break
    return result

def simplified_dfs_expansion(sd: SuccessionDiagram):
    root = sd.root()
    root_space = sd.node_space(root)   # This is after percolation, so it might not be the whole state space.
    root_space_symbolic = symbolic_space(sd.symbolic, root_space)

    if len(root_space) == sd.network.num_vars():
        # A special case where the network has one fixed-point/minimal trap that is
        # resolvable through percolation. Not super common but 
        sd.expanded.add(root)
        return

    root_motifs = trappist(sd.petri_net,  problem="max", ensure_subspace=root_space)
    root_motifs = sorted(root_motifs, key=lambda x: len(x), reverse=True)

    stack: list[tuple[int, list[dict[str, int]], ColoredVertexSet]] = [(root, root_motifs, root_space_symbolic)]

    while len(stack) > 0:
        top = stack[-1]
        node, motifs, remaining = top
        node_space = sd.node_space(node)

        if DEBUG:
            print(f"[{len(stack)}] Process node {node} with {len(motifs)} motifs and {remaining.cardinality()} states.")

        if len(motifs) == 0:
            if DEBUG:
                if remaining.is_empty():                    
                    print(" >> Done.")
                else:
                    print(f" >> Done. Found trap avoidant attractor: {remaining.cardinality()}")
            sd.expanded.add(node)
            stack.pop()
            continue
        
        if remaining.is_empty():            
            for motif in motifs:
                child = sd.ensure_node(node, motif)
                if child not in sd.expanded:
                    # It is possible that this process will re-discover an existing
                    # expanded node. In such case, it would be invalid to mark it as stub.
                    sd.node_is_stub(child, True)
            if DEBUG:
                print(f" >> Done. {len(motifs)} stubs created.")
            sd.expanded.add(node)
            stack.pop()            
            continue

        motif = motifs.pop()
        motif_symbolic = symbolic_space(sd.symbolic, motif)
        child = sd.ensure_node(node, motif)        
        child_space = sd.node_space(child)
        child_space_symbolic = symbolic_space(sd.symbolic, child_space)            

        # If this motif is already covered, we can create a stub node.
        # But we should still remove it from the remaining nodes.
        is_redundant = child_space_symbolic.intersect(remaining).is_empty()
        
        remaining = remaining.minus(motif_symbolic)
        if SYMBOLIC_PRUNING:
            free_variables = [ x for x in sd.network.variables() if sd.network.get_variable_name(x) not in node_space ]
            remaining = trap(sd.symbolic, free_variables, remaining)
        
        stack[-1] = (node, motifs, remaining)        

        if DEBUG:
            print(f" >> Expanding child motif into node {child}. Reduced remaining states: {remaining.cardinality()}.")

        if is_redundant:
            if child not in sd.expanded:            
                sd.node_is_stub(child, True)
                if DEBUG:
                    print(f" >> Created child is a stub.")
            else:
                if DEBUG:
                    print(f" >> Child is redundant and already expanded.")
            continue        

        if child not in sd.expanded:
            if len(child_space) == sd.network.num_vars():
                if DEBUG:
                    print(" >> Child is a fixed-point. Done.")
                sd.expanded.add(child)
            else:
                child_motifs = trappist(sd.petri_net,  problem="max", ensure_subspace=child_space)
                if len(child_motifs) == 0:
                    if DEBUG:
                        print(f" >> Child is a trap ({len(child_space)}/{sd.network.num_vars()} fixed vars). Done.")
                    sd.expanded.add(child)
                else:                
                    root_motifs = sorted(child_motifs, key=lambda x: len(x), reverse=True)
                    stack.append((child, child_motifs, child_space_symbolic))

        



def simplified_bfs_expansion_2(
    sd: SuccessionDiagram,
    node_id: int,
    depth_limit: int | None = None,
    node_limit: int | None = None,
) -> int:
    level_nodes = set([node_id])
    visited = set()
    
    level = 0
    level_remaining = { 0: sd.symbolic.unit_colored_vertices() }

    total_expanded = 0
    while len(level_nodes) > 0:
        print(f" >> Start level with {len(level_nodes)} nodes.")
        next_level = set()

        assert level+1 not in level_remaining
        level_remaining[level+1] = sd.symbolic.empty_colored_vertices()        

        for i, node in enumerate(sorted(level_nodes, key=lambda node: len(sd.node_space(node)), reverse=True)):
            print(f" >> {i+1} / {len(level_nodes)}")
            if node in visited:                
                continue
            visited.add(node)

            node_space = sd.node_space(node)
            node_space_symbolic = symbolic_space(sd.symbolic, node_space)

            # BN variables that are not fixed in this space
            free_variables = [ x for x in sd.network.variables() if sd.network.get_variable_name(x) not in node_space ]

            if node not in sd.expanded:
                if node_space_symbolic.intersect(level_remaining[level]).is_empty():
                    # At this point, all previous nodes have successfully covered the
                    # state space of this node. Hence there's no need to expand it and
                    # it can remain as a stub.
                    sd.node_is_stub(node, True)
                    continue

                assert not sd.node_is_stub(node)

                total_expanded += 1
                sd.expanded.add(node)                

                if len(node_space) == sd.network.num_vars():
                    # This node is a fixed-point. Trappist would just
                    # return this fixed-point again. No need to continue.
                    if DEBUG:
                        print(f"Found fixed-point: {node_space}.")
                    continue

                if DEBUG:
                    print(f"[{node}] Expanding with {len(node_space)} fixed vars.")
            
                sub_spaces = trappist(
                    sd.petri_net, 
                    problem="max", 
                    ensure_subspace=node_space,
                )

                if len(sub_spaces) == 0:
                    if DEBUG:
                        print(f"Found minimum trap space: {node_space}.")
                    continue

                if DEBUG:
                    print(f"Sub-spaces: {len(sub_spaces)}")

                for child in sub_spaces:
                    sd.ensure_node(node, child)
                
            else:
                # A single BFS search should only visit nodes that are not expanded.
                # But if you rerun the search again on a diagram that is partially
                # expanded, expanded nodes can start appearing here. In such case,
                # we want to continue exploring, but there is no need to actually
                # do anything with this node.                                 
                pass

            # At this point, node is either expanded or a stub.
            if not sd.node_is_stub(node):
                remaining = level_remaining[level].minus(node_space_symbolic)
                if SYMBOLIC_PRUNING:
                    remaining = trap(sd.symbolic, free_variables, remaining)
                level_remaining[level] = remaining            

            for s in sd.successors(node):
                if s not in visited and not sd.node_is_stub(s):
                    # Add the successor node into the next level, and add its state
                    # space into the "remaining", uncovered state space of that level.
                    s_space = sd.node_space(s)
                    s_space_symbolic = symbolic_space(sd.symbolic, s_space)
                    level_remaining[level+1] = level_remaining[level+1].union(s_space_symbolic)
                    next_level.add(s)

            if DEBUG:
                expanded_frac = f"{len(sd.expanded)}/{sd.G.number_of_nodes()}"
                expanded_perc = 100 * len(sd.expanded) / sd.G.number_of_nodes()
                print(f"Total expanded: {expanded_frac} ({round(expanded_perc)}%).")

            if (node_limit is not None) and (total_expanded >= node_limit):                    
                return total_expanded

        if level_remaining[level].is_empty():
            print(f"No attractor candidates at level {level}.")
        else:
            print(f"Found attractor candidates at level {level}. Remaining: {level_remaining[level].cardinality()}")

        if (depth_limit is not None) and (level >= depth_limit):
            break

        level += 1
        level_nodes = next_level        

    return total_expanded


def simplified_bfs_expansion(sd: SuccessionDiagram, node_id: int, depth_limit: int | None = None, node_limit: int | None = None) -> int:
    """
    This algorithm is based on the "exhaustive" BFS search implemented within `SuccessionDiagram.expand_node`.
    However, it uses a local criterion to reduce the number of expanded nodes (the excluded nodes are preserved
    as "stub" nodes). 

    The correctness of the method hinges on two observations:
     - First, no attractor or minimal trap space can reside in the "percolated state space".
     That is, states residing in `space \ percolate(space)` are safe to disregard in terms of 
     SD construction, since they don't contain any "significant" features.
     - Second, let C be a child subspace within a specific parent SD node, and S a set of (different)
     fully expanded (i.e. non-stub) child nodes of the same parent. Then, assuming `percolate(C)`
     is the subset of S, expanding C is redundant and C can become a stub. This is because all minimal 
     trap spaces and attractors that reside within C also reside within *some* member of S.

    As such, for every SD node, it is sufficient to expand a subset of its child nodes that "covers"
    the *percolated* state space of the remaining child nodes. Computing a minimal such subset is
    a search problem that I leave as a TODO, since it seems non-trivial and it is not clear how much it would
    actually improve performance. However, a simple greedy algorithm can be implemented in a few lines of Python
    and produces a reasonably small set of "covering child spaces".

    Finally, note that this principle leads to the conclusion that every minimal trap space is reached
    using a path that is equivalent to that in a "full" succession diagram. This is because every normal
    node (not a stub) has exactly the same successors as in the "full" expansion algorithm. Hence, assuming
    a minimal trap space is discovered, it must be discovered using a path that also appears in the full
    succession diagram.
    """

    bfs_queue = [(node_id, depth_limit)]
    visited = set()        

    total_expanded = 0

    while len(bfs_queue) > 0:
        node, depth = bfs_queue.pop(0)

        # Due to BFS, a node is always first visited with its maximal
        # depth_limit, hence there is no need to re-explore nodes, even
        # if they are visited again with a different depth limit.
        if node in visited:
            continue
        visited.add(node)  

        # If the node isn't expanded, we can expand it:
        if node not in sd.expanded:                
            total_expanded += 1

            sd.expanded.add(node)

            # Stubs should be marked as expanded and as such
            # should not appear here.
            assert not sd.G.nodes[node]['stub']
            
            current_space = sd.node_space(node)
            current_space_symbolic = symbolic_space(sd.symbolic, current_space)

            if len(current_space) == sd.network.num_vars():
                # This node is a fixed-point. Trappist would just
                # return this fixed-point again. No need to continue.
                if DEBUG:
                    print(f"Found fixed-point: {current_space}.")
                continue

            if DEBUG:
                print(f"[{node}] Expanding with {len(current_space)} fixed vars.")
            
            sub_spaces = trappist(
                sd.petri_net, 
                problem="max", 
                ensure_subspace=current_space,
            )

            if len(sub_spaces) == 0:
                if DEBUG:
                    print(f"Found minimum trap space: {current_space}.")
                continue

            if DEBUG:
                print(f"Sub-spaces: {len(sub_spaces)}")

            # First, sort sub-spaces such that the largest
            # subspace (based on cardinality) is last.
            sub_spaces = sorted(sub_spaces, key=lambda s: len(s), reverse=True)

            to_expand_symbolic = sd.symbolic.empty_colored_vertices()

            to_expand = []
            to_skip = []

            free_variables = [ x for x in sd.network.variables() if sd.network.get_variable_name(x) not in current_space ]

            while len(sub_spaces) > 0:
                child = sub_spaces.pop()
                child_symbolic = symbolic_space(sd.symbolic, child)
                child_percolated, _ = percolate_space(sd.network, child, strict_percolation=False)
                child_percolated_symolic = symbolic_space(sd.symbolic, child_percolated)

                if child_percolated_symolic.is_subset(to_expand_symbolic):
                    to_skip.append(child)
                else:
                    to_expand.append(child)
                    to_expand_symbolic = to_expand_symbolic.union(child_symbolic)
                    
                    if SYMBOLIC_PRUNING:
                        to_expand_symbolic = reach_bwd(sd.symbolic, free_variables, to_expand_symbolic, current_space_symbolic)

            if current_space_symbolic.is_subset(to_expand_symbolic):
                if DEBUG:
                    print(f"Ruled-out trap avoidant attractors.")
                sd.attr_expanded.add(node)
            else:
                if DEBUG:
                    print(f"Possible trap avoidant attractor here.")

            if DEBUG:
                print(f"{len(to_skip)} sub-spaces skipped, {len(to_expand)} expanded.")

            # TODO: For some reason, model 079 is very slow here. 
            # The only reason I can think of is percolation/PyEDA parsing recursion depth?

            for sub_space in to_skip:
                child_id = sd.ensure_node(node, sub_space)
                sd.node_is_stub(child_id, True)

            for sub_space in to_expand:
                sd.ensure_node(node, sub_space)

            if DEBUG:
                expanded_frac = f"{len(sd.expanded)}/{sd.G.number_of_nodes()}"
                expanded_perc = 100 * len(sd.expanded) / sd.G.number_of_nodes()                                    
                fixed_vars_frac = f"{len(sd.node_space(node))}/{sd.network.num_vars()}"
                print(f"Total expanded: {expanded_frac} ({round(expanded_perc)}%). Fixed vars {fixed_vars_frac} at depth {sd.node_depth(node)}.")

            if (node_limit is not None) and (total_expanded >= node_limit):                    
                return total_expanded
        
        # If the node has sufficient depth limit, we can 
        # explore its successors (as long as they are not visited).

        if (depth is None) or (depth > 0):
            new_depth = None if depth is None else (depth - 1)
            for s in sd.G.successors(node):
                if s not in visited:
                    bfs_queue.append((s, new_depth))

    return total_expanded
