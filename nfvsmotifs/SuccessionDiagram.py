from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Set, Iterator
    from biodivine_aeon import BooleanNetwork

import networkx as nx  # type: ignore

from nfvsmotifs.interaction_graph_utils import feedback_vertex_set
from nfvsmotifs.petri_net_translation import network_to_petrinet
from nfvsmotifs.space_utils import percolate_space, space_unique_key
from nfvsmotifs.trappist_core import trappist

from nfvsmotifs._sd_algorithms.expand_bfs import expand_bfs
from nfvsmotifs._sd_algorithms.expand_dfs import expand_dfs
from nfvsmotifs._sd_algorithms.expand_minimal_spaces import expand_minimal_spaces
from nfvsmotifs._sd_algorithms.compute_attractor_seeds import compute_attractor_seeds
from nfvsmotifs._sd_algorithms.expand_attractor_seeds import expand_attractor_seeds

# Enables helpful "progress" messages.
DEBUG = False


class SuccessionDiagram:
    """
    `SuccessionDiagram` (SD) is a directed acyclic graph representing the structure of trap spaces 
    induced by a particular `BooleanNetwork`. The root of the graph is the whole state space (after
    percolation of constant values). The leaf nodes are individual minimal trap spaces. Each path 
    from the root to a leaf represents a succession of gradually more restrictive trap spaces.

    There are several additional "features" implemented by `SuccessionDiagram` that allow us to
    implement more advanced (and efficient) algorithms:
        - The succession diagram can be expanded lazily. Each node is initially a *stub* node,
        meaning none of its child nodes are known. A stub can be then *expanded* into a full node
        by computing the *stable motifs* for the associated trap space.
        - Each node can be annotated with *attractor seed states*, i.e. states that are known to
        cover all network attractors within that space.        

    Overview of the `SuccessionDiagram` API:
        - Introspection of the whole diagram:
            * `SuccessionDiagram.root()`: the integer ID of the root node.
            * `SuccessionDiagram.depth()`: The depth (hight?) of the current succession diagram.
            * `SuccessionDiagram.node_ids() / .stub_ids() / .expanded_ids()`: Iterates over the 
            corresponding *node IDs* managed by this SD.
            * `SuccessionDiagram.is_subgraph() / .is_isomorphic()`: Compares two succession diagrams.
        - Inspecting individual nodes:
            * `SuccessionDiagram.node_depth(id)`: The depth (length of a maximal path from the 
            root) of the given node.
            * `SuccessionDiagram.node_space(id)`: The trap space associated with the given node.
            * `SuccessionDiagram.node_is_minimal(id)`: Checks if the given node is a minimal
            trap space: i.e. it is expanded and has no successor nodes.
            * `SuccessionDiagram.node_is_expaned(id)`: Check if the given node is expanded, i.e.
            its successor nodes are computed.
            * `SuccessionDiagram.node_successors(id, compute=True/False)`: The list of successor 
            node IDs for the given node (can be computed if not yet known, in which case the
            node becomes expanded).
            * `SuccessionDiagram.node_attractor_seeds(id, compute=True/False)`: The list of 
            "attractor seed states" associated with the given node (if these are not known yet,
            they can be computed).
            * `SuccessionDiagram.edge_stable_motif(id, id)`: Obtain the stable motif which enables
            the edge between the two nodes.
        - There are several "expand procedures" that can explore a larger part of the SD at once,
        typically with a specific goal in mind:
            * `SuccessionDiagram.expand_bfs() / .expand_dfs()`: Expands the whole SD up to a certain
            upper bound on size/depth.
            * `SuccessionDiagram.expand_minimal_traps()`: Expands the SD such that each minimal trap
            space is reachable by at least one path from the root.

    *Other internal implementation notes:* 
    
    *The node IDs are assumed to be a continuous range of integers. If this breaks at any point, 
    please make sure to check the correctness of the new implementation thoroughly.*

    *There is no way to remove or otherwise "contract" the succession diagram: Once a node is expanded, 
    it stays expanded. At the time of writing, there does not appear to be a need for deleting SD 
    nodes and it greatly simplifies reasoning about correctness.*

    """

    def __init__(self, network: BooleanNetwork):
        # Original Boolean network.
        self.network = network
        # A Petri net representation of the original Boolean network.
        self.petri_net = network_to_petrinet(network)
        # Negative feedback vertex set.
        self.nfvs = feedback_vertex_set(
            network, parity="negative"
        )  # find_minimum_NFVS(network)
        # A directed acyclic graph representing the succession diagram.
        self.G = nx.DiGraph()
        # A dictionary used for uniqueness checks on the nodes of the succession diagram.
        # See `SuccessionDiagram.ensure_node` for details.
        self.node_indices: dict[int, int] = {}

        # Create an un-expanded root node.        
        self._ensure_node(None, {})

    def __len__(self) -> int:
        """
        Returns the number of nodes in this `SuccessionDiagram`.
        """
        return self.G.number_of_nodes()

    def root(self) -> int:
        """
        Return the ID of the root node.
        """
        return 0

    def depth(self) -> int:
        """
        Compute the maximal node depth in the diagram.

        Depth is counted from zero (root has depth zero).
        """
        d = 0
        for node in self.G.nodes():  # type: ignore[reportUnknownVariableType] # noqa
            d = max(d, self.node_depth(int(node)))  # type: ignore[reportUnknownArgumentType] # noqa
        return d

    def node_ids(self) -> Iterator[int]:
        """
        Iterator over all available node IDs.
        """
        for i in range(len(self)):
            yield i

    def stub_ids(self) -> Iterator[int]:
        """
        Iterator over all node IDs that are currently not expanded.
        """
        for i in range(len(self)):
            if not self.node_is_expanded(i):
                yield i

    def expanded_ids(self) -> Iterator[int]:
        """
        Iterator over all node IDs that are currently expanded.
        """
        for i in range(len(self)):
            if self.node_is_expanded(i):
                yield i

    def minimal_trap_spaces(self) -> list[int]:
        """
        List of node IDs that represent the minimal trap spaces within this succession diagram.

        Note that stub nodes do not count as minimal!
        """
        return [ i for i in self.expanded_ids() if self.node_is_minimal(i) ]

    def find_node(self, node_space: dict[str, int]) -> int | None:
        """
        Return the ID of the node matching the provided `node_space`, or `None` if no such
        node exists in this succession diagram.
        """
        try:
            key = space_unique_key(node_space, self.network)
            if key in self.node_indices:
                i = self.node_indices[key]
                # This assertion could be violated if a user gives a node space that is not based
                # on the same network as this succession diagram.
                assert node_space == self.node_space(i)
                return i
            else:
                return None
        except RuntimeError:
            # If the user gives us a space that uses variables not used by this network, 
            # we should get an error that we can catch and report that no such space exists here.
            return None                    

    def is_subgraph(self, other: SuccessionDiagram) -> bool:
        """
        Returns `True` if this succession diagram is a subgraph of the `other` succession diagram.
        
        Note that this function works even for diagrams based on different Boolean
        networks, as long as both succession diagrams only depend on the same subset of 
        network variables.

        WARNING: This does not take into account the stable motifs on individual edges. Just the
        subspaces associated with nodes and the presence of edges between nodes.
        """        
        # Every stub node is reachable through an expanded node and 
        # thus will be checked by the following code.
        for i in self.expanded_ids():
            other_i = other.find_node(self.node_space(i))
            if other_i is None:
                return False
            my_successors = self.node_successors(i)
            other_successors = []
            if other.node_is_expanded(other_i):                
                other_successors = other.node_successors(other_i)

            for my_s in my_successors:
                other_s = other.find_node(self.node_space(my_s))
                if other_s not in other_successors:
                    return False
        return True

    def is_isomorphic(self, other: SuccessionDiagram) -> bool:
        """
        Returns `True` if the two succession diagrams are isomorphic.

        Note that this function works even for diagrams based on different Boolean
        networks, as long as both succession diagrams only depend on the same subset of 
        network variables.

        WARNING: This does not take into account the stable motifs on individual edges. Just the
        subspaces associated with nodes and the presence of edges between nodes.
        """
        return self.is_subgraph(other) and other.is_subgraph(self)

    def node_depth(self, node_id: int) -> int:
        """
        Get the depth associated with the provided `node_id`. The depth is counted
        as the longest path from the root node to the given node.
        """        
        return self.G.nodes[node_id]["depth"]  # type: ignore[reportUnknownVariableType] # noqa

    def node_space(self, node_id: int) -> dict[str, int]:
        """
        Get the sub-space associated with the provided `node_id`.

        Note that this is the space *after* percolation. Hence it can hold that
        `|node_space(child)| < |node_space(parent)| + |stable_motif(parent, child)|`.
        """
        return self.G.nodes[node_id]["space"]  # type: ignore[reportUnknownVariableType] # noqa

    def node_is_expanded(self, node_id: int) -> bool:
        """
        True if the successors of the given node are already computed.
        """
        return self.G.nodes[node_id]["expanded"]

    def node_is_minimal(self, node_id: int) -> bool:
        """
        True if the node is expanded and it has no successors, i.e. it is a minimal trap space.
        """
        is_leaf: bool = self.G.out_degree(node_id) == 0  # type: ignore[reportUnknownMemberType, reportUnknownVariableType] # noqa
        return is_leaf and self.G.nodes[node_id]["expanded"]  # type: ignore[reportUnknownVariableType] # noqa

    def node_successors(self, node_id: int, compute: bool = False) -> list[int]:
        """
        Return the successor nodes for the given `node_id`. If the node is already expanded, known results
        are simply returned. If the node is not expanded, but `compute` is set to `True`, then the node
        is expanded and the newly computed results are returned. If the node is not expanded and `compute` 
        is set to `False`, the method raises a `KeyError` exception.

        The default behaviour intentionally does not compute successors to prevent "accidental complexity".

        WARNING: We do not guarantee the order of the returned nodes. If you need the successors in a 
        repeatable order, you should sort the list manually.

        Also note that if the given `node_id` already has associated attractor data but is not expanded,
        this data will be deleted as it is no longer up to date.
        """
        node = self.G.nodes[node_id]
        if not node["expanded"] and not compute:
            raise KeyError(f"Node {node_id} is not expanded.")
        if not node["expanded"]:
            self._expand_one_node(node_id)
        
        return list(self.G.successors(node_id))

    def node_attractor_seeds(self, node_id: int, compute: bool = False) -> list[dict[str, int]]:
        """
        Return the list of attractor seed states corresponding to the given `node_id`. Similar to
        `node_successors`, the method either computes the data if unknown, or throws an exception,
        depending on the `compute` flag.

        Note that you can compute attractor seeds for stub nodes, but (a) these attractors are not
        guaranteed to be unique (i.e. you can "discover" the same attractor in multiple stub nodes,
        if the stub nodes intersect), and (b) this data is erased if the stub node is expanded 
        later on.
        """
        node = self.G.nodes[node_id]

        attractors = node["attractors"]
        
        if attractors is None and not compute:
            raise KeyError(f"Attractor data not computed for node {node_id}.")
        
        if attractors is None:
            attractors = compute_attractor_seeds(self, node_id)
            node["attractors"] = attractors            

        return attractors 

    def edge_stable_motif(self, parent_id: int, child_id: int) -> dict[str, int]:
        """
        Return the *stable motif* associated with the specified parent-child edge.

        This corresponds to the maximal trap space within the `parent_id` node that, after percolation,
        yields the `child_id` node.
        """
        return self.G.edges[parent_id, child_id]["motif"]  # type: ignore[reportUnknownVariableType] # noqa

    def expand_bfs(
        self, 
        node_id: int | None = None, 
        bfs_level_limit: int | None = None, 
        size_limit: int | None = None,
    ) -> bool:
        """
        Explore the succession diagram in a BFS manner.
            - If `node_id` is given, initiate BFS from this node. Otherwise use root.
            - If `bfs_level_limit` is given, this is the last "level" (distance from the initial node) 
            of nodes that should be expanded (any subsequent child nodes are left unexplored).
            - If `size_limit` is given, the procedure stops once `SuccessionDiagram` exceeds the given size.

        With default settings, the method will explore the whole succession diagram without any restrictions.

        The method returns `True` if the whole exploration was completed and `False` if it was terminated
        early based on one of the aforementioned conditions.

        Note that the procedure also explores nodes that are already expanded. I.e. if all nodes at levels
        0,1,2 are expanded, but there are stub nodes on level 3, the procedure will still discover and 
        expand these stub nodes (assuming sufficient level and size limit).

        Also note that the `size_limit` is only a soft limit: for each node, we always have to create all
        child nodes when expanding it. Hence the procedure can only check the condition in between
        expanding new nodes.
        """
        return expand_bfs(self, node_id, bfs_level_limit, size_limit)

    def expand_dfs(
        self,
        node_id: int | None = None,
        dfs_stack_limit: int | None = None,
        size_limit: int | None = None,
    ) -> bool:
        """
        Similar to `expand_bfs`, but uses DFS instead of BFS.

        The only major difference is the `dfs_stack_limit` which restricts the size of the DFS stack.
        Nodes that would appear "deeper" in the stack than this limit are left unexpanded. Note that 
        this stack size is technically *some* form of distance from the initial node, but not necessarily
        the minimal distance.
        """
        return expand_dfs(self, node_id, dfs_stack_limit, size_limit)

    def expand_minimal_spaces(self):
        """
        Expands the succession diagram in a way that guarantees every minimal trap space to be reachable 
        from the root node, but otherwise (greedily) avoids unnecesary expansion of nodes whenever possible.

        The algorithm is loosely based on `expand_bfs` implementation, but on each BFS level only expands the
        first node that still contains some minimal trap space not covered by a previously expanded node 
        at that level.

        The resulting succession diagram construction is deterministic, but can vary if some nodes are already 
        expanded initially. In such case, the procedure still tries to avoid expanding unnecessary nodes,
        which means existing expanded nodes can be prioritised over the "canonical" ones.
        """
        expand_minimal_spaces(self)

    def expand_attractor_seeds(self):
        """
        Expands the succession diagram such that for every asynchronous attractor, there is
        at least one expanded trap space which is the minimal trap space containing this attractor.
        In other words, the procedure expands the succession diagram as little as possible, but 
        ensures that every attractor is "as easy to identify" as possible.

        After this procedure, it is sufficient to search for attractors in expanded nodes.
        Note that this method does not perform exact attractor identification. It is possible
        that some nodes are expanded spuriously and the succession diagram is thus larger 
        than necessary.
        """
        expand_attractor_seeds(self)

    def _update_node_depth(self, node_id: int, parent_id: int):
        """
        An internal method that updates the depth of a node based on a specific parent node. 
        This assumes that there is an edge from `parent` to `node_id`.

        Note that the depth can only increase.
        """
        assert self.G.edges[parent_id, node_id] is not None
        parent_depth = self.G.nodes[parent_id]["depth"]
        current_depth = self.G.nodes[node_id]["depth"]
        self.G.nodes[node_id]["depth"] = max(current_depth, parent_depth + 1)

    def _expand_one_node(self, node_id: int):
        """
        An internal method that expands a single node of the succession diagram.

        This entails computing the maximal trap spaces within the node (stable motifs)
        and creating a node for the result (if it does not exist yet).

        If the node is already expanded, the method does nothing.

        If there are already some attractor data for this node (stub nodes can have 
        associated attractor data), this data is erased.
        """
        node = self.G.nodes[node_id]
        if node["expanded"]:
            return

        node["expanded"] = True        
        node["attractors"] = None
        current_space = node["space"]

        if DEBUG:
            print(f"[{node_id}] Expanding: {len(self.node_space(node_id))} fixed vars.")

        if len(current_space) == self.network.num_vars():
            # This node is a fixed-point. Trappist would just
            # return this fixed-point again. No need to continue.
            if DEBUG:
                print(f"[{node_id}] Found fixed-point: {current_space}.")
            return

        sub_spaces = trappist(
            self.petri_net,
            problem="max",
            ensure_subspace=current_space,
        )

        # Sort the spaces based on a unique key in case trappist is not always sorted deterministically.
        sub_spaces = sorted(sub_spaces, key=lambda space: space_unique_key(space, self.network))

        if len(sub_spaces) == 0:
            if DEBUG:
                print(f"[{node_id}] Found minimum trap space: {current_space}.")
            return

        if DEBUG:
            print(f"[{node_id}] Found sub-spaces: {len(sub_spaces)}")

        for sub_space in sub_spaces:
            child_id = self._ensure_node(node_id, sub_space)  # type: ignore[reportUnusedVariable] # noqa

            if DEBUG:
                print(f"[{node_id}] Created edge into node {child_id}.")

    def _ensure_node(self, parent_id: int | None, stable_motif: dict[str, int]) -> int:
        """
        Internal method that ensures the provided node is present in this succession diagram as a 
        child of the given `parent_id`.

        The `stable_motif` is an "initial" trap space that is then percolated to compute the actual 
        fixed variables for this node. The method also updates the depth of the child node if necessary.

        If the `parent_id` is not given, no edge is created and depth is considered to be zero
        (i.e. the node is the root).
        """

        fixed_vars, _ = percolate_space(
            self.network, stable_motif, strict_percolation=False
        )

        key = space_unique_key(fixed_vars, self.network)

        child_id = None
        if key not in self.node_indices:
            child_id = self.G.number_of_nodes()
            self.G.add_node(
                child_id,
                id = child_id,    # Just in case we ever need it within the "node data" dictionary. 
                space = fixed_vars,
                depth = 0,
                expanded = False,
                attractors = None,
            )
            self.node_indices[key] = child_id
        else:
            child_id = self.node_indices[key]

        assert child_id is not None
        
        if parent_id is not None:
            # TODO: It seems that there are some networks where the same child can be reached
            # through multiple stable motifs. Not sure how to approach these... but this is 
            # probably good enough for now.
            self.G.add_edge(parent_id, child_id, motif=stable_motif)
            self._update_node_depth(child_id, parent_id)

        return child_id
