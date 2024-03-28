from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal, cast

if TYPE_CHECKING:
    from typing import Iterator

import networkx as nx  # type: ignore
from biodivine_aeon import AsynchronousGraph, BooleanNetwork, VariableId

from balm._sd_algorithms.compute_attractor_seeds import compute_attractor_seeds
from balm._sd_algorithms.expand_attractor_seeds import expand_attractor_seeds
from balm._sd_algorithms.expand_bfs import expand_bfs
from balm._sd_algorithms.expand_dfs import expand_dfs
from balm._sd_algorithms.expand_minimal_spaces import expand_minimal_spaces
from balm._sd_algorithms.expand_source_SCCs import expand_source_SCCs
from balm._sd_algorithms.expand_to_target import expand_to_target
from balm.interaction_graph_utils import (
    cleanup_network,
    feedback_vertex_set,
    source_SCCs,
)
from balm.petri_net_translation import (
    extract_source_variables,
    network_to_petrinet,
    restrict_petrinet_to_subspace,
)
from balm.space_utils import percolate_network, percolate_space, space_unique_key
from balm.trappist_core import trappist
from balm.types import BooleanSpace, NodeData, SuccessionDiagramState

DEBUG: bool = False
"""
If `True`, print debug and progress messages.
"""


class SuccessionDiagram:
    """
    Succession diagram of a Boolean network.

    This encodes relationships between minimal trap spaces and can be used for
    attractor detection and control. Note that the succession diagram is
    expanded lazily, so it is not built until `build` or a similar method is
    called.


    Examples
    --------
    >>> import balm
    >>> sd = balm.SuccessionDiagram.from_rules(\"""
    ...     A, B
    ...     B, A & C
    ...     C, !A | B
    ... \""")
    >>> print(sd.summary()) # not built yet!
    Succession Diagram with 1 nodes and depth 0.
    State order: A, B, C
    <BLANKLINE>
    Attractors in diagram:
    <BLANKLINE>
    >>> sd.build()
    >>> print(sd.summary()) # now it's built
    Succession Diagram with 3 nodes and depth 1.
    State order: A, B, C
    <BLANKLINE>
    Attractors in diagram:
    <BLANKLINE>
    minimal trap space 001
    ...................001
    <BLANKLINE>
    minimal trap space 111
    ...................111
    <BLANKLINE>
    """

    NFVS_NODE_THRESHOLD: int = 2_000
    """
    If the number of nodes in the (redunced) network is greater than this
    threshold, we find a feedback vertex set (without consideration of cycle
    sign) for reachability preprocessing. There is a trade-off between the speed
    gains from a smaller node set to consider and the cost of determining which
    FVS nodes only intersect negative cycles to find an NFVS subset. Typically,
    for smaller networks, the trade-off is in favor of computing a smaller NFVS.
    """

    __slots__ = (
        "network",
        "symbolic",
        "petri_net",
        "nfvs",
        "dag",
        "node_indices",
    )

    def __init__(self, network: BooleanNetwork):
        # Original Boolean network.
        self.network: BooleanNetwork = cleanup_network(network)
        """
        The Boolean network as an AEON network.
        """

        self.symbolic: AsynchronousGraph = AsynchronousGraph(self.network)
        """
        The symbolic representation of the network.
        """

        self.petri_net: nx.DiGraph = network_to_petrinet(network)
        """
        The Petri net representation of the network (see :mod:`petri_net_translation<balm.petri_net_translation>`).
        """

        if DEBUG:
            print(
                f"Generated global Petri net with {len(self.petri_net.nodes)} nodes and {len(self.petri_net.edges)} edges."
            )

        # Initially, we don't need the NFVS for anything.
        self.nfvs: list[str] | None = None
        """
        The negative feedback vertex set used for attractor detection.
        """

        self.dag: nx.DiGraph = nx.DiGraph()
        """
        The directed acyclic graph (DAG) representation of the succession
        diagram structure.
        """
        # A dictionary used for uniqueness checks on the nodes of the succession
        # diagram. See `SuccessionDiagram.ensure_node` for details.
        self.node_indices: dict[int, int] = {}
        """
        A dictionary mapping subspace keys to their positions in the succession
        diagram (see :func:`balm.space_utils.space_unique_key`).
        """

        # Create an un-expanded root node.
        self._ensure_node(None, {})

    def __getstate__(self) -> SuccessionDiagramState:
        return {
            "network_rules": self.network.to_aeon(),
            "petri_net": self.petri_net,
            "nfvs": self.nfvs,
            "dag": self.dag,
            "node_indices": self.node_indices,
        }

    def __setstate__(self, state: SuccessionDiagramState):
        # In theory, the network should be cleaned-up at this point, but just in case...
        self.network = cleanup_network(BooleanNetwork.from_aeon(state["network_rules"]))
        self.symbolic = AsynchronousGraph(self.network)
        self.petri_net = state["petri_net"]
        self.nfvs = state["nfvs"]
        self.dag = state["dag"]
        self.node_indices = state["node_indices"]

    def __len__(self) -> int:
        """
        Returns the number of nodes in this `SuccessionDiagram`.
        """
        return self.dag.number_of_nodes()

    @staticmethod
    def from_rules(
        rules: str,
        format: Literal["bnet", "aeon", "sbml"] = "bnet",
    ) -> SuccessionDiagram:
        """
        Generate a succession diagram from the given string.

        Parameters
        ----------
        rules : str
            The string representation of the network rules.
        format : Literal['bnet', 'aeon', 'sbml']
            The format of the string. One of `"bnet"`, `"aeon"`, or `"sbml"`.
            Defaults to `"bnet"`.

        Returns
        -------
        SuccessionDiagram
            The generated succession diagram. Initially unexpanded.
        """

        if format == "bnet":
            return SuccessionDiagram(BooleanNetwork.from_bnet(rules))
        elif format == "aeon":
            return SuccessionDiagram(BooleanNetwork.from_aeon(rules))
        elif format == "sbml":
            return SuccessionDiagram(BooleanNetwork.from_sbml(rules))
        else:
            raise ValueError(f"Unknown format: {format}")

    @staticmethod
    def from_file(path: str) -> SuccessionDiagram:
        """
        Read a `BooleanNetwork` from the given file path. The format is automatically inferred from
        the file extension.
        """
        return SuccessionDiagram(BooleanNetwork.from_file(path))

    def expanded_attractor_seeds(self) -> dict[int, list[BooleanSpace]]:
        """
        Attractor seeds for each expanded node.

        An attractor seed is a state belonging to an attractor (and thus a state
        from which the entire attractor is reachable, by definition).

        If called before the `SuccessionDiagram` is fully built, this will not
        be a complete accounting of attractor seed states.

        Returns
        -------
        dict[int,list[BooleanSpace]]
            Each key is the id of an expanded succession diagram node, whose
            corresponding value is a list of attractor seeds for that node. Note
            that one succession diagram node can have multiple attractors. Ordering
            of the lists in the returned dictionary is not guaranteed.

        Example
        -------
        >>> import balm
        >>> sd = balm.SuccessionDiagram.from_rules(\"""
        ...     A, B
        ...     B, A & C
        ...     C, !A | B
        ... \""")
        >>> sd.build()
        >>> eas = sd.expanded_attractor_seeds()
        >>> for id, atts in sorted(eas.items()):
        ...     for x in atts:
        ...         print(f"{id}: {dict(sorted(x.items()))}")
        1: {'A': 0, 'B': 0, 'C': 1}
        2: {'A': 1, 'B': 1, 'C': 1}
        """
        res: dict[int, list[BooleanSpace]] = {}
        for id in self.expanded_ids():
            atts = self.node_attractor_seeds(id)
            if not atts:  # no attractors for this node
                continue
            res[id] = atts

        return res

    def summary(self) -> str:
        """
        Return a summary of the succession diagram as a string.
        """
        var_ordering = sorted(
            [self.network.get_variable_name(v) for v in self.network.variables()]
        )
        report_string = (
            f"Succession Diagram with {len(self)} nodes and depth {self.depth()}.\n"
            f"State order: {', '.join(var_ordering)}\n\n"
            "Attractors in diagram:\n\n"
        )
        for node in self.node_ids():
            try:
                attrs = self.node_attractor_seeds(node, compute=False)
            except KeyError:
                continue

            if not attrs:
                continue

            space = self.node_data(node)["space"]

            if self.node_is_minimal(node):
                space_str_prefix = "minimal trap space "
            else:
                space_str_prefix = "motif avoidance in "
            space_str = ""
            for var in var_ordering:
                if var in space:
                    space_str += str(space[var])
                else:
                    space_str += "*"
            report_string += f"{space_str_prefix}{space_str}\n"
            for attr in attrs:
                attr_str = "".join(str(v) for _, v in sorted(attr.items()))
                report_string += "." * len(space_str_prefix) + f"{attr_str}\n"
            report_string += "\n"
        # remove final extra newline and return
        return report_string[:-1]

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
        for node in cast(set[int], self.dag.nodes()):
            d = max(d, self.node_data(int(node))["depth"])
        return d

    def node_ids(self) -> Iterator[int]:
        """
        Iterator over all available node IDs.
        """
        for i in range(len(self)):
            yield i

    def stub_ids(self) -> Iterator[int]:
        """
        Iterator over all node IDs that are currently *not* expanded.
        """
        for i in range(len(self)):
            if not self.node_data(i)["expanded"]:
                yield i

    def expanded_ids(self) -> Iterator[int]:
        """
        Iterator over all node IDs that are currently expanded.
        """
        for i in range(len(self)):
            if self.node_data(i)["expanded"]:
                yield i

    def minimal_trap_spaces(self) -> list[int]:
        """
        List of node IDs that represent the minimal trap spaces within this
        succession diagram.

        Note that stub nodes do not count as minimal!
        """
        return [i for i in self.expanded_ids() if self.node_is_minimal(i)]

    def find_node(self, node_space: BooleanSpace) -> int | None:
        """
        Return the ID of the node matching the provided `node_space`, or `None`
        if no such node exists in this succession diagram.

        Parameters
        ----------
        node_space: BooleanSpace
            The space of the node to find.

        Returns
        -------
        int | None
            The ID of the node matching the provided `node_space`, or `None`
            if no such node exists in this succession diagram.
        """
        try:
            key = space_unique_key(node_space, self.network)  # throws IndexError
            if key in self.node_indices:
                return self.node_indices[key]
            else:
                return None
        except IndexError:
            # If `space_unique_key` finds a variable that does not exist in this
            # `SuccessionDiagram`, it throws an `IndexError`. This can happen
            # for example if we are comparing two succession diagrams based on
            # completely different networks.
            return None

    def is_subgraph(self, other: SuccessionDiagram) -> bool:
        """
        Returns `True` if this succession diagram is a subgraph of the `other`
        succession diagram.

        Note that this function works even for diagrams based on different
        Boolean networks, as long as both succession diagrams only depend on the
        same subset of network variables.

        WARNING: This does not take into account the stable motifs on individual
        edges. Just the subspaces associated with nodes and the presence of
        edges between nodes.

        Parameters
        ----------
        other: SuccessionDiagram
            The other succession diagram.

        Returns
        -------
        bool
            `True` if this succession diagram is a subgraph of the `other`
            succession diagram.
        """
        # Every stub node is reachable through an expanded node and
        # thus will be checked by the following code.
        for i in self.expanded_ids():
            other_i = other.find_node(self.node_data(i)["space"])
            if other_i is None:
                return False
            my_successors = self.node_successors(i)
            other_successors = []
            if other.node_data(other_i)["expanded"]:
                other_successors = other.node_successors(other_i)

            for my_s in my_successors:
                other_s = other.find_node(self.node_data(my_s)["space"])
                if other_s not in other_successors:
                    return False
        return True

    def is_isomorphic(self, other: SuccessionDiagram) -> bool:
        """
        Returns `True` if the two succession diagrams are isomorphic.

        Note that this function works even for diagrams based on different
        Boolean networks, as long as both succession diagrams only depend on the
        same subset of network variables.

        WARNING: This does not take into account the stable motifs on individual
        edges. Just the subspaces associated with nodes and the presence of
        edges between nodes.

        Parameters
        ----------
        other: SuccessionDiagram
            The other succession diagram.

        Returns
        -------
        bool
            `True` if the two succession diagrams are isomorphic.
        """
        return self.is_subgraph(other) and other.is_subgraph(self)

    def node_data(self, node_id: int) -> NodeData:
        """
        Get the data associated with the provided `node_id`.

        Returns a `NodeData` object with the following attributes:

        - `depth`: The depth of the node.
        - `attractors`: The attractors of the node (if computed).
        - `petri_net`: The Petri net representation of the node (if computed).
        - `space`: The sub-space of the node.
        - `expanded`: Whether the node is expanded or not.

        See :class:`balm.types.NodeData` for more information.

        Parameters
        ----------
        node_id: int
            The ID of the node.

        Returns
        -------
        NodeData
            The data associated with the provided `node_id`. Note that at
            runtime, this object is an untyped dictionary.
        """
        return cast(NodeData, self.dag.nodes[node_id])

    def node_is_minimal(self, node_id: int) -> bool:
        """
        True if the node represents a minimal trap space.

        Parameters
        ----------
        node_id: int
            The ID of the node.

        Returns
        -------
        bool
            `True` if the node is expanded and it has no successors, i.e. it is a
            minimal trap space.
        """
        is_leaf: bool = self.dag.out_degree(node_id) == 0  # type: ignore
        return is_leaf and self.dag.nodes[node_id]["expanded"]  # type: ignore

    def node_successors(self, node_id: int, compute: bool = False) -> list[int]:
        """
        Return the successor nodes for the given `node_id`.

        If the node is already expanded, known results are simply returned. If
        the node is not expanded, but `compute` is set to `True`, then the node
        is expanded and the newly computed results are returned. If the node is
        not expanded and `compute` is set to `False`, the method raises a
        `KeyError` exception.

        The default behaviour intentionally does not compute successors to
        prevent "accidental complexity".

        WARNING: We do not guarantee the order of the returned nodes. If you
        need the successors in a repeatable order, you should sort the list
        manually.

        Also note that if the given `node_id` already has associated attractor
        data but is not expanded, this data will be deleted as it is no longer
        up to date.

        Parameters
        ----------
        node_id: int
            The ID of the node.
        compute: bool
            Whether to compute the successors if they are not already known.

        Returns
        -------
        list[int]
            The list of successor node ids.
        """
        node = cast(dict[str, Any], self.dag.nodes[node_id])
        if not node["expanded"] and not compute:
            raise KeyError(f"Node {node_id} is not expanded.")
        if not node["expanded"]:
            self._expand_one_node(node_id)

        return list(self.dag.successors(node_id))  # type: ignore

    def node_attractor_seeds(
        self, node_id: int, compute: bool = False
    ) -> list[BooleanSpace]:
        """
        Return the list of attractor seed states for the given `node_id`.

        Similar to :meth:`node_successors`, the method either computes the
        data if unknown, or throws an exception, depending on the `compute`
        flag.

        Note that you can compute attractor seeds for stub nodes, but (a) these
        attractors are not guaranteed to be unique (i.e. you can "discover" the
        same attractor in multiple stub nodes, if the stub nodes intersect), and
        (b) this data is erased if the stub node is expanded later on.

        Parameters
        ----------
        node_id: int
            The ID of the node.
        compute: bool
            Whether to compute the attractor seeds if they are not already known.

        Returns
        -------
        list[BooleanSpace]
            The list of attractor seed states.
        """
        node = cast(NodeData, self.dag.nodes[node_id])

        attractors = node["attractors"]

        if attractors is None and not compute:
            raise KeyError(f"Attractor data not computed for node {node_id}.")

        if attractors is None:
            attractors = compute_attractor_seeds(self, node_id)
            node["attractors"] = attractors

        return attractors

    def node_nfvs(self, node_id: int) -> list[str]:
        """
        Approximate minimum negative feedback vertex set on the subspace for the given node.

        See :func:`balm.interaction_graph_utils.feedback_vertex_set` for further details.

        Parameters
        ----------
        node_id: int
            The ID of the node.

        Returns
        -------
        list[str]
            The negative feedback vertex set, as a list of node names.
        """

        # TODO: Right now, we are just returning the gloval FVS. But in the future,
        # we probably want to compute a separate FVS for each node, because it could
        # be much smaller if a lot of the nodes are fixed.
        assert node_id in self.dag.nodes

        if self.nfvs is None:
            if self.network.variable_count() < SuccessionDiagram.NFVS_NODE_THRESHOLD:
                # Computing the *negative* variant of the FVS is surprisingly costly.
                # Hence it mostly makes sense for the smaller networks only.
                self.nfvs = feedback_vertex_set(self.network, parity="negative")
            else:
                self.nfvs = feedback_vertex_set(self.network)

        return self.nfvs

    def edge_stable_motif(
        self, parent_id: int, child_id: int, reduced: bool = False
    ) -> BooleanSpace:
        """
        Return the stable motif for the given parent-child edge.

        If `reduced` is set to `False` (default), the unpercolated stable
        motif trap space corresponding to the child node is returned; this
        includes the nodes that are fixed in the percolated trap space of the
        parent node. If `reduced` is set to `True`, the nodes fixed in the
        parent are removed (and thus the reduced stable motif is not a trap
        space of the original network, but is a maximal trap space in the
        network reduced by the parent node).

        Parameters
        ----------
        parent_id: int
            The ID of the parent node.
        child_id: int
            The ID of the child node.
        reduced: bool
            Whether to return the reduced stable motif.

        Returns
        -------
        BooleanSpace
            The stable motif (maximal trap space) represented by the edge.
        """

        if reduced:
            return cast(
                BooleanSpace,
                {
                    k: v
                    for k, v in self.dag.edges[parent_id, child_id]["motif"].items()  # type: ignore
                    if k not in self.node_data(parent_id)["space"]
                },
            )
        else:
            return cast(BooleanSpace, self.dag.edges[parent_id, child_id]["motif"])

    def component_subdiagrams(
        self,
        node_id: int | None = None,
    ) -> Iterator[SuccessionDiagram]:
        """
        Return unexpanded subdiagrams for the source SCCs in a node subspace.

        The subnetwork on which the subdiagram is defined is defined by the
        variables in `component_variables`, which is a list of variable names.
        The `component_variables` must be backward-closed, meaning there is no
        variable outside this list that regulates any variable in the
        subnetwork. Note that this is not explicitly checked in this function.

        Also note that the symbolic encoding of the new network is not
        compatible with the encoding of the original network, because the
        underlying networks have different sets of variables.

        Parameters
        ----------
        node_id : int | None
            The ID of a succession diagram node that will define a subspace on
            which the subnetworks should be considered. By default, the root node
            is used.

        Returns
        -------
        Iterator[SuccessionDiagram]
            An iterator over unexpanded succession diagrams of the subnetwork.
        """

        if node_id is None:
            node_id = self.root()

        reference_bn = percolate_network(
            self.network,
            self.node_data(node_id)["space"],
            remove_constants=True,
        )

        source_scc_list = source_SCCs(reference_bn)

        for component_variables in source_scc_list:
            new_bn = BooleanNetwork(component_variables)

            # Build a mapping between the old and new network variables.
            id_map: dict[VariableId, VariableId] = {}
            for var in component_variables:
                old_id = reference_bn.find_variable(var)
                assert old_id is not None
                new_id = new_bn.find_variable(var)
                assert new_id is not None
                id_map[old_id] = new_id

            # Copy regulations that are in the source component.
            for reg in reference_bn.regulations():
                if reg["source"] in id_map and reg["target"] in id_map:
                    new_bn.add_regulation(
                        {
                            "source": reference_bn.get_variable_name(reg["source"]),
                            "target": reference_bn.get_variable_name(reg["target"]),
                            "essential": reg["essential"],
                            "sign": reg["sign"],
                        }
                    )

            # Copy update functions from the source component after translating them
            # to the new IDs.
            for var_id in id_map.keys():
                old_function = reference_bn.get_update_function(var_id)
                assert old_function is not None
                new_function = old_function.rename_all(new_bn, variables=id_map)
                new_bn.set_update_function(id_map[var_id], new_function)

            yield SuccessionDiagram(new_bn)

    def build(self):
        """
        Expand the succession diagram and search for attractors using default methods.
        """
        self.expand_scc()
        for node_id in self.node_ids():
            self.node_attractor_seeds(node_id, compute=True)

    def expand_scc(self, find_motif_avoidant_attractors: bool = True) -> bool:
        """
        Expand the succession diagram using the source SCC method.
        """
        return expand_source_SCCs(self, check_maa=find_motif_avoidant_attractors)

    def expand_bfs(
        self,
        node_id: int | None = None,
        bfs_level_limit: int | None = None,
        size_limit: int | None = None,
    ) -> bool:
        """
        Explore the succession diagram in a BFS manner.

        If `node_id` is given, initiate BFS from this node. Otherwise use root.
        If `bfs_level_limit` is given, this is the last "level" (distance from
        the initial node) of nodes that should be expanded (any subsequent child
        nodes are left unexplored).

        If `size_limit` is given, the procedure stops once `SuccessionDiagram`
        exceeds the given size.

        With default settings, the method will explore the whole succession
        diagram without any restrictions.

        The method returns `True` if the whole exploration was completed and
        `False` if it was terminated early based on one of the aforementioned
        conditions.

        Note that the procedure also explores nodes that are already expanded.
        I.e. if all nodes at levels 0,1,2 are expanded, but there are stub nodes
        on level 3, the procedure will still discover and expand these stub
        nodes (assuming sufficient level and size limit).

        Also note that the `size_limit` is only a soft limit: for each node, we
        always have to create all child nodes when expanding it. Hence the
        procedure can only check the condition in between expanding new nodes.
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

        The only major difference is the `dfs_stack_limit` which restricts the
        size of the DFS stack. Nodes that would appear "deeper" in the stack
        than this limit are left unexpanded. Note that this stack size is
        technically *some* form of distance from the initial node, but not
        necessarily the minimal distance.
        """
        return expand_dfs(self, node_id, dfs_stack_limit, size_limit)

    def expand_minimal_spaces(self, size_limit: int | None = None) -> bool:
        """
        Expands the succession diagram in a way that guarantees every minimal
        trap space to be reachable from the root node, but otherwise (greedily)
        avoids unnecesary expansion of nodes whenever possible.

        The algorithm is loosely based on `expand_bfs` implementation, but on
        each BFS level only expands the first node that still contains some
        minimal trap space not covered by a previously expanded node at that
        level.

        The resulting succession diagram construction is deterministic, but can
        vary if some nodes are already expanded initially. In such case, the
        procedure still tries to avoid expanding unnecessary nodes, which means
        existing expanded nodes can be prioritised over the "canonical" ones.
        """
        return expand_minimal_spaces(self, size_limit)

    def expand_attractor_seeds(self, size_limit: int | None = None) -> bool:
        """
        Expands the succession diagram such that for every asynchronous
        attractor, there is at least one expanded trap space which is the
        minimal trap space containing this attractor. In other words, the
        procedure expands the succession diagram as little as possible, but
        ensures that every attractor is "as easy to identify" as possible.

        After this procedure, it is sufficient to search for attractors in
        expanded nodes. Note that this method does not perform exact attractor
        identification. It is possible that some nodes are expanded spuriously
        and the succession diagram is thus larger than necessary.
        """
        return expand_attractor_seeds(self, size_limit)

    def expand_to_target(
        self, target: BooleanSpace, size_limit: int | None = None
    ) -> bool:
        """
        Expands the succession diagram using BFS in such a way that only nodes
        which intersect `target` but are not fully contained in it are expanded.

        This is used for example in control, as it ensures that all possible
        branches of the succession diagram relevant for a particular "target
        subspace" are expanded as much as necessary, but not more.
        """
        return expand_to_target(self, target, size_limit)

    def _update_node_depth(self, node_id: int, parent_id: int):
        """
        An internal method that updates the depth of a node based on a specific
        parent node. This assumes that there is an edge from `parent` to
        `node_id`.

        Note that the depth can only increase.
        """
        assert self.dag.edges[parent_id, node_id] is not None
        parent_depth = cast(int, self.dag.nodes[parent_id]["depth"])
        current_depth = cast(int, self.dag.nodes[node_id]["depth"])
        self.dag.nodes[node_id]["depth"] = max(current_depth, parent_depth + 1)

    def _update_node_petri_net(self, node_id: int, parent_id: int | None):
        """
        Try to compute a restricted Petri net for this node using the Petri net
        stored in the node's parent.

        If the parent Petri net does not exist or the parent is not specified, the
        restriction is performed on the "global" Petri net object, which is usually
        slower but should yield the same result.

        No Petri net is computed if the specified node is a fixed-point, since
        such Petri net is always empty.
        """

        node_space = self.node_data(node_id)["space"]

        if len(node_space) == self.network.variable_count():
            # If fixed point, no need to compute.
            return

        if parent_id is not None:
            parent_pn = self.node_data(parent_id)["petri_net"]
            if parent_pn is None:
                pn = self.petri_net
            else:
                pn = parent_pn
        else:
            pn = self.petri_net

        restricted_pn = restrict_petrinet_to_subspace(pn, node_space)

        if DEBUG:
            print(
                f"[{node_id}] Generated Petri net restriction with {len(restricted_pn.nodes)} nodes and {len(restricted_pn.edges)} edges."
            )

        self.dag.nodes[node_id]["petri_net"] = restricted_pn

    def _clear_node_petri_net(self, node_id: int):
        """
        Remove the computed Petri net for this node.

        This is typically done once the node is expanded, since we know the Petri net won't
        be needed anymore.
        """
        del self.dag.nodes[node_id]["petri_net"]

    def _expand_one_node(self, node_id: int):
        """
        An internal method that expands a single node of the succession diagram.

        This entails computing the maximal trap spaces within the node (stable
        motifs) and creating a node for the result (if it does not exist yet).

        If the node is already expanded, the method does nothing.

        If there are already some attractor data for this node (stub nodes can
        have associated attractor data), this data is erased.
        """
        node = cast(dict[str, Any], self.dag.nodes[node_id])
        if node["expanded"]:
            return

        node["expanded"] = True
        node["attractors"] = None
        current_space = node["space"]

        if DEBUG:
            print(
                f"[{node_id}] Expanding: {len(self.node_data(node_id)['space'])} fixed vars."
            )

        if len(current_space) == self.network.variable_count():
            # This node is a fixed-point. Trappist would just
            # return this fixed-point again. No need to continue.
            if DEBUG:
                print(f"[{node_id}] Found fixed-point: {current_space}.")
            return

        # We use the non-propagated Petri net for backwards-compatibility reasons here.
        # The SD created from the restricted Petri net is technically correct, but can
        # propagate some of the input values further and yields a smaller SD.
        source_nodes = []
        if node_id == self.root():
            source_nodes = extract_source_variables(self.petri_net)

        sub_spaces: list[BooleanSpace]
        pn = self.node_data(node_id)["petri_net"]
        if pn is not None:
            # We have a pre-propagated PN for this sub-space, hence we can use
            # that to compute the trap spaces.
            partial_sub_spaces = trappist(
                pn, problem="max", optimize_source_variables=source_nodes
            )
            sub_spaces = [(s | current_space) for s in partial_sub_spaces]
        else:
            # If we (for whatever reason) don't have the pre-propagated PN,
            # we can still use the "global" PN and let trappist deal with the restriction.
            sub_spaces = trappist(
                self.petri_net,
                problem="max",
                ensure_subspace=current_space,
                optimize_source_variables=source_nodes,
            )

        # Sort the spaces based on a unique key in case trappist is not always
        # sorted deterministically.
        sub_spaces = sorted(
            sub_spaces, key=lambda space: space_unique_key(space, self.network)
        )

        if len(sub_spaces) == 0:
            if DEBUG:
                print(f"[{node_id}] Found minimum trap space: {current_space}.")
            return

        if DEBUG:
            print(f"[{node_id}] Found sub-spaces: {len(sub_spaces)}")

        for sub_space in sub_spaces:
            child_id = self._ensure_node(node_id, sub_space)

            if DEBUG:
                print(f"[{node_id}] Created edge into node {child_id}.")

        # The node is now fully expanded, hence we won't need the restricted PN anymore.
        self._clear_node_petri_net(node_id)

    def _ensure_node(self, parent_id: int | None, stable_motif: BooleanSpace) -> int:
        """
        Internal method that ensures the provided node is present in this
        succession diagram as a child of the given `parent_id`.

        The `stable_motif` is an "initial" trap space that is then percolated to
        compute the actual fixed variables for this node. The method also
        updates the depth of the child node if necessary.

        If the `parent_id` is not given, no edge is created and depth is
        considered to be zero (i.e. the node is the root).
        """

        fixed_vars = percolate_space(self.symbolic, stable_motif)

        key = space_unique_key(fixed_vars, self.network)

        child_id = None
        if key not in self.node_indices:
            child_id = self.dag.number_of_nodes()

            # Note: this must match the fields of the `NodeData` class
            self.dag.add_node(  # type: ignore
                child_id,
                space=fixed_vars,
                depth=0,
                expanded=False,
                attractors=None,
            )
            self.node_indices[key] = child_id
        else:
            child_id = self.node_indices[key]

        assert child_id is not None

        if parent_id is not None:
            # TODO: It seems that there are some networks where the same child
            # can be reached through multiple stable motifs. Not sure how to
            # approach these... but this is probably good enough for now.
            self.dag.add_edge(parent_id, child_id, motif=stable_motif)  # type: ignore
            self._update_node_depth(child_id, parent_id)

        self._update_node_petri_net(child_id, parent_id)

        return child_id
