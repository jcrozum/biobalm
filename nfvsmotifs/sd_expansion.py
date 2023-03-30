from __future__ import annotations

from nfvsmotifs.SuccessionDiagram import SuccessionDiagram
from nfvsmotifs.space_utils import percolate_space
from nfvsmotifs.trappist_core import trappist

DEBUG = False
SYMBOLIC_PRUNING = True

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
                if DEBUG and result.symbolic_size() > 100_000:                    
                    print(f"BDD node count: {int(result.symbolic_size() / 10)}; Reached {result.cardinality()}/{universe.cardinality()}")
                keep_working = True

                # TODO: This limit should by somehow proportional to the size of the
                # network... something like 100_000 * variable_count seems reasonable?
                # Furthermore, we might want to have a limit on the total number of iterations.
                if result.symbolic_size() > 1_000_000:
                    return result
                break
    return result


def symbolic_subspace(
    stg: SymbolicAsyncGraph,
    subspace: dict[str, int]
) -> ColoredVertexSet:
    transformed = { x:bool(y) for x,y in subspace.items() }
    return stg.fix_subspace(transformed)

def simplified_bfs_expansion(sd: SuccessionDiagram, node_id: int, depth_limit: int | None = None, node_limit: int | None = None):
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
            current_space_symbolic = symbolic_subspace(sd.symbolic, current_space)

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
                child_symbolic = symbolic_subspace(sd.symbolic, child)
                child_percolated, _ = percolate_space(sd.network, child, strict_percolation=False)
                child_percolated_symolic = symbolic_subspace(sd.symbolic, child_percolated)

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
