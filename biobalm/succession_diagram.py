from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal, cast

if TYPE_CHECKING:
    from typing import Iterator

import copy
import networkx as nx  # type: ignore
from biodivine_aeon import AsynchronousGraph, BooleanNetwork, VertexSet

# Attractor detection algorithms.
from biobalm._sd_attractors.attractor_candidates import compute_attractor_candidates
from biobalm._sd_attractors.attractor_symbolic import (
    compute_attractors_symbolic,
    symbolic_attractor_fallback,
)

# SD expansion algorithms/heuristics.
from biobalm._sd_algorithms.expand_attractor_seeds import expand_attractor_seeds
from biobalm._sd_algorithms.expand_bfs import expand_bfs
from biobalm._sd_algorithms.expand_dfs import expand_dfs
from biobalm._sd_algorithms.expand_minimal_spaces import expand_minimal_spaces
from biobalm._sd_algorithms.expand_source_SCCs import expand_source_SCCs
from biobalm._sd_algorithms.expand_source_blocks import expand_source_blocks
from biobalm._sd_algorithms.expand_to_target import expand_to_target

from biobalm.interaction_graph_utils import (
    cleanup_network,
    feedback_vertex_set,
    source_SCCs,
)
from biobalm.petri_net_translation import (
    extract_source_variables,
    network_to_petrinet,
    restrict_petrinet_to_subspace,
)
from biobalm.space_utils import (
    percolate_network,
    percolate_space,
    space_unique_key,
    is_subspace,
)
from biobalm.trappist_core import trappist
from biobalm.types import (
    BooleanSpace,
    NodeData,
    SuccessionDiagramState,
    SuccessionDiagramConfiguration,
)


class SuccessionDiagram:
    """
    Succession diagram of a Boolean network.

    This encodes relationships between minimal trap spaces and can be used for
    attractor detection and control. Note that the succession diagram is
    expanded lazily, so it is not built until `build` or a similar method is
    called.


    Examples
    --------
    >>> import biobalm
    >>> sd = biobalm.SuccessionDiagram.from_rules(\"""
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

    __slots__ = (
        "network",
        "symbolic",
        "petri_net",
        "nfvs",
        "dag",
        "node_indices",
        "config",
    )

    def __init__(
        self,
        network: BooleanNetwork,
        config: SuccessionDiagramConfiguration | None = None,
    ):
        if config is None:
            config = SuccessionDiagram.default_config()
        self.config = config

        # Original Boolean network.
        self.network: BooleanNetwork = cleanup_network(network)
        """
        The Boolean network represented as a `biodivine_aeon.BooleanNetwork` object.
        """

        self.symbolic: AsynchronousGraph = AsynchronousGraph(self.network)
        """
        The symbolic representation of the network using `biodivine_aeon.AsynchronousGraph`.
        """

        self.petri_net: nx.DiGraph = network_to_petrinet(network)
        """
        The Petri net representation of the network (see :mod:`petri_net_translation<biobalm.petri_net_translation>`).
        """

        if self.config["debug"]:
            print(
                f"Generated global Petri net with {len(self.petri_net.nodes)} nodes and {len(self.petri_net.edges)} edges."
            )

        self.nfvs: list[str] | None = None
        """
        The negative feedback vertex set of `SuccessionDiagram.network`, or `None` if not computed yet.
        """

        self.dag: nx.DiGraph = nx.DiGraph()
        """
        The directed acyclic graph (DAG) representation of the succession
        diagram structure.
        """

        self.node_indices: dict[int, int] = {}
        """
        A dictionary mapping subspace keys to their positions in the succession
        diagram (see :func:`biobalm.space_utils.space_unique_key`).
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
            "config": self.config,
        }

    def __setstate__(self, state: SuccessionDiagramState):
        # In theory, the network should be cleaned-up at this point, but just in case...
        self.network = cleanup_network(BooleanNetwork.from_aeon(state["network_rules"]))
        self.symbolic = AsynchronousGraph(self.network)
        self.petri_net = state["petri_net"]
        self.nfvs = state["nfvs"]
        self.dag = state["dag"]
        self.node_indices = state["node_indices"]
        self.config = state["config"]

    def __len__(self) -> int:
        """
        Returns the number of nodes in this `SuccessionDiagram`.
        """
        return self.dag.number_of_nodes()

    @staticmethod
    def default_config() -> SuccessionDiagramConfiguration:
        return {
            "debug": False,
            "max_motifs_per_node": 100_000,
            "nfvs_size_threshold": 2_000,
            "pint_goal_size_limit": 8_192,
            "attractor_candidates_limit": 100_000,
            "retained_set_optimization_threshold": 1_000,
            "minimum_simulation_budget": 1_000,
        }

    @staticmethod
    def from_rules(
        rules: str,
        format: Literal["bnet", "aeon", "sbml"] = "bnet",
        config: SuccessionDiagramConfiguration | None = None,
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
        config : SuccessionDiagramConfiguration | None
            An optional configuration object with internal settings
            and default values.

        Returns
        -------
        SuccessionDiagram
            The generated succession diagram. Initially unexpanded.
        """

        if format == "bnet":
            return SuccessionDiagram(BooleanNetwork.from_bnet(rules), config)
        elif format == "aeon":
            return SuccessionDiagram(BooleanNetwork.from_aeon(rules), config)
        elif format == "sbml":
            return SuccessionDiagram(BooleanNetwork.from_sbml(rules), config)
        else:
            raise ValueError(f"Unknown format: {format}")

    @staticmethod
    def from_file(
        path: str, config: SuccessionDiagramConfiguration | None = None
    ) -> SuccessionDiagram:
        """
        Read a `BooleanNetwork` from the given file path. The format is automatically inferred from
        the file extension.

        Optionally, you can also supply a configuration object to customize the
        resulting succession diagram.
        """
        return SuccessionDiagram(BooleanNetwork.from_file(path), config)

    def expanded_attractor_candidates(self) -> dict[int, list[BooleanSpace]]:
        """
        Attractor candidates for each expanded node. The candidate list is
        computed for nodes where it is not known yet.

        Collectively, for every attractor in every expanded node, this returns
        at least one candidate state for said attractor. However, for complex
        attractors, there can be more than one candidate state, and there can
        be candidates that are not contained in any attractor.

        *If called before the `SuccessionDiagram` is fully built, this will not
        be a complete accounting of attractors, since any node that isn't expanded
        is not included in the result.*

        See also:
         - :meth:`expanded_attractor_seeds<SuccessionDiagram.expanded_attractor_seeds>`
         - :meth:`expanded_attractor_seeds<SuccessionDiagram.expanded_attractor_sets>`
         - :meth:`expanded_attractor_seeds<SuccessionDiagram.node_attractor_candidates>`
         - :meth:`expanded_attractor_seeds<SuccessionDiagram.node_attractor_seeds>`
         - :meth:`expanded_attractor_seeds<SuccessionDiagram.node_attractor_sets>`

        Returns
        -------
        dict[int,list[BooleanSpace]]
            Each key is the id of an expanded succession diagram node, whose
            corresponding value is a list of attractor candidates for that node. Note
            that one succession diagram node can have multiple attractors. Ordering
            of the lists in the returned dictionary is not guaranteed.

        Example
        -------
        >>> import biobalm
        >>> sd = biobalm.SuccessionDiagram.from_rules(\"""
        ...     A, B
        ...     B, A & C
        ...     C, !A | B
        ... \""")
        >>> sd.build()
        >>> eas = sd.expanded_attractor_candidates()
        >>> for id, atts in sorted(eas.items()):
        ...     for x in atts:
        ...         print(f"{id}: {dict(sorted(x.items()))}")
        1: {'A': 0, 'B': 0, 'C': 1}
        2: {'A': 1, 'B': 1, 'C': 1}
        """
        res: dict[int, list[BooleanSpace]] = {}
        for id in self.expanded_ids():
            atts = self.node_attractor_candidates(id, compute=True)
            if not atts:  # no attractors for this node
                continue
            res[id] = atts

        return res

    def expanded_attractor_seeds(self) -> dict[int, list[BooleanSpace]]:
        """
        Attractor seeds for each expanded node. The seed list is
        computed for nodes where it is not known yet.

        Collectively, for every attractor in every expanded node, this returns
        exactly one seed state for said attractor.

        *If called before the `SuccessionDiagram` is fully built, this will not
        be a complete accounting of attractors, since any node that isn't expanded
        is not included in the result.*

        See also:
         - :meth:`expanded_attractor_seeds<SuccessionDiagram.expanded_attractor_candidates>`
         - :meth:`expanded_attractor_seeds<SuccessionDiagram.expanded_attractor_sets>`
         - :meth:`expanded_attractor_seeds<SuccessionDiagram.node_attractor_candidates>`
         - :meth:`expanded_attractor_seeds<SuccessionDiagram.node_attractor_seeds>`
         - :meth:`expanded_attractor_seeds<SuccessionDiagram.node_attractor_sets>`

        Returns
        -------
        dict[int,list[BooleanSpace]]
            Each key is the id of an expanded succession diagram node, whose
            corresponding value is a list of attractor seeds for that node. Note
            that one succession diagram node can have multiple attractors. Ordering
            of the lists in the returned dictionary is not guaranteed.

        Example
        -------

        *Note that for this simple network, attractor candidates and attractor
        seeds are the same states.*

        >>> import biobalm
        >>> sd = biobalm.SuccessionDiagram.from_rules(\"""
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
            atts = self.node_attractor_seeds(id, compute=True)
            if not atts:  # no attractors for this node
                continue
            res[id] = atts

        return res

    def expanded_attractor_sets(self) -> dict[int, list[VertexSet]]:
        """
        Attractor sets for each expanded node. The sets are
        computed for nodes where they are not known yet.

        These sets represent the complete collection of attractor states in each
        expanded node. For fixed-point attractors, this is effectively equivalent
        to the attractor seeds. For complex attractors, this set contains the
        full attractor. Hence, it is harder to compute but can facilitate richer
        post-processing and analysis.

        *If called before the `SuccessionDiagram` is fully built, this will not
        be a complete accounting of attractors, since any node that isn't expanded
        is not included in the result.*

        See also:
         - :meth:`expanded_attractor_seeds<SuccessionDiagram.expanded_attractor_candidates>`
         - :meth:`expanded_attractor_seeds<SuccessionDiagram.expanded_attractor_seeds>`
         - :meth:`expanded_attractor_seeds<SuccessionDiagram.node_attractor_candidates>`
         - :meth:`expanded_attractor_seeds<SuccessionDiagram.node_attractor_seeds>`
         - :meth:`expanded_attractor_seeds<SuccessionDiagram.node_attractor_sets>`

        Returns
        -------
        dict[int,list[biodivine_aeon.VertexSet]]
            Each key is the id of an expanded succession diagram node, whose
            corresponding value is a list of attractor sets for that node. Note
            that one succession diagram node can have multiple attractors. Ordering
            of the lists in the returned dictionary is not guaranteed.

        Example
        -------
        >>> import biobalm
        >>> sd = biobalm.SuccessionDiagram.from_rules(\"""
        ...     A, B
        ...     B, A & C
        ...     C, !A | B
        ... \""")
        >>> sd.build()
        >>> eas = sd.expanded_attractor_sets()
        >>> for id, atts in sorted(eas.items()):
        ...     for x in atts:
        ...         print(f"{id}: {x}")
        1: VertexSet(cardinality=1, symbolic_size=5)
        2: VertexSet(cardinality=1, symbolic_size=5)
        """
        res: dict[int, list[VertexSet]] = {}
        for id in self.expanded_ids():
            atts = self.node_attractor_sets(id, compute=True)
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
        - `space`: The sub-space of the node.
        - `expanded`: Whether the node is expanded or not.
        - `percolated_network`: [`None` if not computed] The percolation of
        `SuccessionDiagram.network` with respect to the node's `space`.
        - `percolated_petri_net`: [`None` if not computed] The percolation of
        `SuccessionDiagram.petri_net` with respect to the node's `space`.
        - `percolated_nfvs`: [`None` if not computed] The NFVS of `percolated_network`.
        - `attractor_candidates`: [`None` if not computed] A collection of states
        that collectively cover every attractor of this node.
        - `attractor_seeds`: [`None` if not computed] A collection of states
        one-to-one corresponding to the attractors of this node.
        - `attractor_sets`: [`None` if not computed] A complete collection of
        attractors of this node.

        See :class:`biobalm.types.NodeData` for more information.

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

    def reclaim_node_data(self):
        """
        Removes non-essential data from the `NodeData` dictionary of each node.

        This method can be used to reduce the memory footprint of the succession
        diagram, especially before serialization. However, note that this can
        also slow down subsequent computations if the erased data needs
        to be re-computed.

        The method removes the `percolated_network`, `percolated_petri_net`,
        and the `percolated_nfvs`. Furthermore, if `attractor_seeds` are known,
        it erases the `attractor_candidates`, since seeds can be used for the
        same tasks.
        """

        for node_id in self.node_ids():
            data = self.node_data(node_id)
            data["percolated_network"] = None
            data["percolated_petri_net"] = None
            data["percolated_nfvs"] = None
            if data["attractor_seeds"] is not None:
                data["attractor_candidates"] = None

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
        node = self.node_data(node_id)
        if not node["expanded"] and not compute:
            raise KeyError(f"Node {node_id} is not expanded.")
        if not node["expanded"]:
            self._expand_one_node(node_id)

        return list(self.dag.successors(node_id))  # type: ignore

    def node_attractor_candidates(
        self,
        node_id: int,
        compute: bool = False,
        greedy_asp_minification: bool = True,
        simulation_minification: bool = True,
        pint_minification: bool = False,
    ) -> list[BooleanSpace]:
        """
        Return the list of attractor candidate states for the given `node_id`.

        If attractor candidates are not computed but seeds are, returns
        attractor seeds, as these are also valid as candidate states, but
        even more precise.

        Similar to :meth:`node_successors`, the method either computes the
        data if unknown, or throws an exception, depending on the `compute`
        flag. If `compute` is set to `True`, additional flags can be used
        to adjust the candidate identification process (see *Parameters*).

        Note that you can compute attractor candidates for nodes that are not expanded,
        but (a) multiple unexpanded nodes can contain the same attractor, and hence also
        the same/similar candidates (i.e. you can "discover" the same attractor in multiple
        unexpanded nodes, if the nodes intersect), and (b) this data is erased if the
        node is later expanded.

        Parameters
        ----------
        node_id: int
            The ID of the node.
        compute: bool
            Whether to compute the attractor candidates if they are not already known.
        greedy_asp_minification: bool
            Indicate that the initial candidate set should be first greedily minified
            through repeated ASP queries. [Default: True]
        simulation_minification: bool
            Indicate that the candidate set should be refined through stochastic
            simulation. [Default: True]
        pint_minification: bool
            Indicate that the candidate set should be refined through reachability
            analysis using `pint`. Only enable this option if you actually have `pint`
            installed (it is an optional dependency). [Default: False]

        Returns
        -------
        list[BooleanSpace]
            The list of attractor candidate states.
        """
        node = self.node_data(node_id)

        candidates = node["attractor_candidates"]

        # If candidates are already cleared, but seeds are known, we can
        # just use seeds.
        if candidates is None and node["attractor_seeds"] is not None:
            return node["attractor_seeds"]

        if candidates is None and not compute:
            raise KeyError(f"Attractor candidates not computed for node {node_id}.")

        if candidates is None:
            candidates = compute_attractor_candidates(
                self,
                node_id,
                greedy_asp_minification,
                simulation_minification,
                pint_minification,
            )
            node["attractor_candidates"] = candidates

            # If the computed candidates are actually valid as seeds, just
            # propagate the value so that it doesn't need to be computed later.
            node_is_pseudo_minimal = (not node["expanded"]) or self.node_is_minimal(
                node_id
            )
            if len(candidates) == 0 or (
                node_is_pseudo_minimal and len(candidates) == 1
            ):
                node["attractor_seeds"] = candidates

        return candidates

    def node_attractor_seeds(
        self,
        node_id: int,
        compute: bool = False,
        symbolic_fallback: bool = False,
    ) -> list[BooleanSpace]:
        """
        Return the list of attractor seed states for the given `node_id`.

        Similar to :meth:`node_successors`, the method either computes the
        data if unknown, or throws an exception, depending on the `compute`
        flag.

        Note that the same considerations regarding attractors in unexpanded
        nodes apply as for :meth:`node_attractor_candidates`.

        Parameters
        ----------
        node_id: int
            The ID of the node.
        compute: bool
            Whether to compute the attractor seeds if they are not already known.
        symbolic_fallback: bool
            If active, the method will attempt to compute the attractor seeds fully
            symbolically if the default NFVS-based method fails. However, note that
            the program can become unresponsive if the symbolic encoding of the
            result grows to be too large. [Default: False]

        Returns
        -------
        list[BooleanSpace]
            The list of attractor seed states.
        """
        node = self.node_data(node_id)

        seeds = node["attractor_seeds"]

        if seeds is None and not compute:
            raise KeyError(f"Attractor seeds not computed for node {node_id}.")

        try:
            if seeds is None:
                candidates = self.node_attractor_candidates(node_id, compute=True)
                # Typically, this should be done when computing the candidates, but just in case
                # something illegal happended... if we can show that the current candidate set
                # is optimal, we just keep it and don't compute the attractors symbolically.
                node_is_pseudo_minimal = (not node["expanded"]) or self.node_is_minimal(
                    node_id
                )
                if len(candidates) == 0 or (
                    node_is_pseudo_minimal and len(candidates) == 1
                ):
                    node["attractor_seeds"] = candidates
                    seeds = candidates
                else:
                    result = compute_attractors_symbolic(
                        self, node_id, candidate_states=candidates, seeds_only=True
                    )
                    node["attractor_seeds"] = result[0]
                    # At this point, attractor_sets could be `None`, but that
                    # is valid, as long as we actually compute them later when
                    # they are needed.
                    node["attractor_sets"] = result[1]
                    seeds = result[0]

                # Release memory once attractor seeds are known. We might need these
                # for attractor set computation later, but only if the seeds are not empty.
                if len(seeds) == 0:
                    node["percolated_network"] = None
                    node["percolated_nfvs"] = None
                    node["percolated_petri_net"] = None
        except RuntimeError as e:
            if not symbolic_fallback:
                raise e

            # The NFVS method failed, likely because the candidate set was too large.
            # We can still try to fix this and compute the attractors symbolically.
            (seeds, sets) = symbolic_attractor_fallback(self, node_id)

            node["attractor_seeds"] = seeds
            node["attractor_sets"] = sets

            if len(seeds) == 0:
                node["percolated_network"] = None
                node["percolated_nfvs"] = None
                node["percolated_petri_net"] = None

        return seeds

    def node_attractor_sets(
        self,
        node_id: int,
        compute: bool = False,
    ) -> list[VertexSet]:
        """
        Return the list of attractor sets for the given `node_id`.

        Similar to :meth:`node_successors`, the method either computes the
        data if unknown, or throws an exception, depending on the `compute`
        flag.

        Note that the same considerations regarding attractors in unexpanded
        nodes apply as for :meth:`node_attractor_candidates`.

        Parameters
        ----------
        node_id: int
            The ID of the node.
        compute: bool
            Whether to compute the attractor sets if they are not already known.

        Returns
        -------
        list[biodivine_aeon.VertexSet]
            The list of attractor sets.
        """
        node = self.node_data(node_id)

        sets = node["attractor_sets"]

        if sets is None and not compute:
            raise KeyError(f"Attractor sets not computed for node {node_id}.")

        if sets is None:
            seeds = self.node_attractor_seeds(node_id, compute=True)
            result: tuple[list[BooleanSpace], list[VertexSet] | None] = ([], [])
            if len(seeds) > 0:
                result = compute_attractors_symbolic(
                    self, node_id, candidate_states=seeds
                )
            assert result[1] is not None
            node["attractor_sets"] = result[1]
            sets = result[1]

        return sets

    def node_percolated_nfvs(self, node_id: int, compute: bool = False) -> list[str]:
        """
        Approximate minimum negative feedback vertex set on the Boolean network
        percolated to the node's sub-space.

        Similar to :meth:`node_successors`, the method either computes the
        data if unknown, or throws an exception, depending on the `compute`
        flag.

        See :func:`biobalm.interaction_graph_utils.feedback_vertex_set` for
        further details.

        Parameters
        ----------
        node_id: int
            The ID of the node.
        compute: bool
            Whether to compute the node NFVS if it is not already known.

        Returns
        -------
        list[str]
            The negative feedback vertex set, as a list of node names.
        """

        assert node_id in self.dag.nodes

        node = self.node_data(node_id)

        if node["percolated_nfvs"] is None and not compute:
            raise KeyError(f"NFVS not computed for node {node_id}.")

        if node["percolated_nfvs"] is None:
            percolated_network = self.node_percolated_network(node_id, compute)
            percolated_size = percolated_network.variable_count()
            if percolated_size < self.config["nfvs_size_threshold"]:
                # Computing the *negative* variant of the FVS is surprisingly costly.
                # Hence it mostly makes sense for the smaller networks only.
                nfvs = feedback_vertex_set(percolated_network, parity="negative")
            else:
                nfvs = feedback_vertex_set(percolated_network)
            node["percolated_nfvs"] = nfvs
        else:
            nfvs = node["percolated_nfvs"]

        return nfvs

    def node_percolated_network(
        self, node_id: int, compute: bool = False
    ) -> BooleanNetwork:
        """
        The Boolean network percolated to the node's sub-space, with
        constant variables removed.

        Similar to :meth:`node_successors`, the method either computes the
        data if unknown, or throws an exception, depending on the `compute`
        flag.

        Parameters
        ----------
        node_id: int
            The ID of the node.
        compute: bool
            Whether to compute the node NFVS if it is not already known.

        Returns
        -------
        biodivine_aeon.BooleanNetwork
            The percolated Boolean network.
        """

        assert node_id in self.dag.nodes

        node = self.node_data(node_id)
        network = node["percolated_network"]

        node_space = node["space"]

        if len(node_space) == self.network.variable_count():
            # If fixed point, no need to compute, the network is always empty.
            return BooleanNetwork()

        if network is None and not compute:
            raise KeyError(f"Percolated network not computed for node {node_id}.")

        if network is None:
            network = percolate_network(
                self.network, node_space, self.symbolic, remove_constants=True
            )
            if self.config["debug"]:
                print(
                    f"[{node_id}] Computed percolated network with {network.variable_count()} variables (vs {self.network.variable_count()})."
                )
            node["percolated_network"] = network

        return network

    def node_percolated_petri_net(
        self,
        node_id: int,
        compute: bool = False,
        parent_id: int | None = None,
    ) -> nx.DiGraph:
        """
        The Petri net representation of the Boolean network percolated to the
        node's sub-space (with constant variables removed).

        Similar to :meth:`node_successors`, the method either computes the
        data if unknown, or throws an exception, depending on the `compute`
        flag.

        Parameters
        ----------
        node_id: int
            The ID of the node.
        compute: bool
            Whether to compute the node NFVS if it is not already known.
        parent_id: int | None
            If given, the percolation process starts with the Petri net of the given
            parent node (if computed). If parent is not given, the percolation starts
            with `SuccessionDiagram.petri_net`, which can be slower but yields the
            same result.

        Returns
        -------
        networkx.DiGraph
            The percolated Boolean network.
        """

        assert node_id in self.dag.nodes

        node = self.node_data(node_id)
        percolated_pn = node["percolated_petri_net"]

        node_space = node["space"]

        if len(node_space) == self.network.variable_count():
            # If fixed point, the result is always empty.
            return nx.DiGraph()

        if percolated_pn is None and not compute:
            raise KeyError(f"Percolated network not computed for node {node_id}.")

        if parent_id is None:
            # If no parent node is directly provided, we can try to use a cached one.
            parent_id = node["parent_node"]

        if percolated_pn is None:
            base_pn = self.petri_net
            percolate_space = node_space
            if parent_id is not None:
                parent_pn = self.node_data(parent_id)["percolated_petri_net"]
                if parent_pn is not None:
                    base_pn = parent_pn
                    percolate_space = node_space

            percolated_pn = restrict_petrinet_to_subspace(base_pn, percolate_space)

            if self.config["debug"]:
                print(
                    f"[{node_id}] Generated Petri net restriction with {len(percolated_pn.nodes)} nodes and {len(percolated_pn.edges)} edges."
                )

            node["percolated_petri_net"] = percolated_pn

        return percolated_pn

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

    def component_subdiagram(
        self,
        component_variables: list[str],
        node_id: int | None = None,
    ) -> SuccessionDiagram:
        """
        Return an *unexpanded* `SuccessionDiagram` that is restricted to
        a subnetwork induced by the provided `component_variables`.  Furthermore,
        If `node_id` is given, the subnetwork is first percolated to the
        subspace of the specified node.

        The `component_variables` must be backward-closed in the considered network
        (i.e. either the full network, or the percolated network if `node_id` is given),
        meaning there is no variable outside this list that regulates any variable in the
        subnetwork. If this is not satisfied, the function will fail while
        creating the subnetwork.

        Also note that the symbolic encoding of the new network is not
        compatible with the encoding of the original network, because the
        underlying networks have different sets of variables.

        Parameters
        ----------
        component_variables : list[str]
            Names of variables which induce the subnetwork of the resulting
            succession diagram.
        node_id : int | None
            The ID of a succession diagram node that will define a subspace
            to which the subnetwork is percolated. If not given, the full
            network is considered.

        Returns
        -------
        SuccessionDiagram
            An unexpanded succession diagram of the subnetwork.
        """

        network = self.network
        if node_id is not None:
            network = self.node_percolated_network(node_id, compute=True)

        to_remove = [
            v for v in network.variable_names() if v not in component_variables
        ]
        component_bn = network.drop(to_remove)
        config_copy: SuccessionDiagramConfiguration = copy.copy(self.config)
        return SuccessionDiagram(component_bn, config_copy)

    def source_scc_subdiagrams(
        self,
        node_id: int | None = None,
    ) -> Iterator[SuccessionDiagram]:
        """
        Return unexpanded subdiagrams for the source SCCs in a node subspace.

        Note that the symbolic encoding of the new network is not
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

        reference_bn = self.node_percolated_network(node_id, compute=True)
        source_scc_list = source_SCCs(reference_bn)

        for component_variables in source_scc_list:
            yield self.component_subdiagram(component_variables, node_id)

    def build(self):
        """
        Expand the succession diagram and search for attractors using default methods.
        """
        self.expand_block()
        for node_id in self.node_ids():
            self.node_attractor_seeds(node_id, compute=True)

    def expand_scc(self, find_motif_avoidant_attractors: bool = True) -> bool:
        """
        Expand the succession diagram using the source SCC method.
        """
        return expand_source_SCCs(self, check_maa=find_motif_avoidant_attractors)

    def expand_block(
        self,
        find_motif_avoidant_attractors: bool = True,
        size_limit: int | None = None,
        optimize_source_nodes: bool = True,
        exact_attractor_detection: bool = False,
    ) -> bool:
        """
        Expand the succession diagram using the source block method.

        There is a minor difference in behavior depending on `find_motif_avoidant_attractors`.
        If set to `False`, the expansion only expands one "source block" for each node,
        without checking any attractor properties. If set to `True`, the expansion might
        expand some nodes fully to uncover nodes that precisely cover motif
        avoidant attractors. As a byproduct, if set to `True` and no motif avoidant attractors
        are detected for some node, this is result is saved and the attractors don't
        need to be recomputed later.

        By default, the method also detects any source nodes and directly expands these
        into trap spaces where all source nodes are fixed. This has no correctness impact on
        attractor search and always produces a smaller succession diagram, but if you need to
        obtain a succession diagram where this does not happen (e.g. for testing), you can turn
        this off using `optimize_source_nodes`.

        If `exact_attractor_detection` is selected, the method will use attractor seeds instead
        of attractor candidates to check for motif-avoidant attractors. This means that the attractor
        detection can take much longer. In particular, it will not expand the SD if it cannot detect
        all attractors, meaning it can get "stuck". However, assuming you want to run full attractor
        detection anyway, this might save you some time, as you won't need to re-run  the failed
        candidate detection. (This is only relevant for models where the attractors are very
        complex and the candidate state detection can fail with default settings).
        """
        return expand_source_blocks(
            self,
            find_motif_avoidant_attractors,
            size_limit=size_limit,
            optimize_source_nodes=optimize_source_nodes,
            check_maa_exact=exact_attractor_detection,
        )

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

    def expand_minimal_spaces(
        self,
        node_id: int | None = None,
        size_limit: int | None = None,
        skip_ignored: bool = False,
    ) -> bool:
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

        Optionally, you can start the expansion from a specific node that is not
        the root using `node_id`, or restrict the size of the succession diagram
        with `size_limit`.

        Returns `True` if the expansion procedure terminated without exceeding
        the size limit.

        If `skip_ignored` is set, any nodes that are not expanded by this procedure
        are "skipped" with edges redirected to the corresponding minimal trap spaces
        (see also `SuccessionDiagram.skip_remaining`). This is usually faster than
        using `skip_remaining` directly, since the minimal trap spaces are only
        computed once. However, note that if used with `size_limit`, nodes that
        are not processed when `size_limit` is reached remain unprocessed
        (i.e. it is not guaranteed that all nodes are either fully expanded or skipped).
        """
        return expand_minimal_spaces(self, node_id, size_limit, skip_ignored)

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

    def skip_to_minimal(self, node_id: int) -> bool:
        """
        Skip the expansion of this node (see also `NodeData.skipped`) and add extra
        edges that connect it directly to its minimal trap spaces. Returns `False`
        if the node is already expanded.

        Note that this method is relatively inefficient when applied to multiple
        nodes repeatedly, as it has to recompute the minimal trap spaces for each node.
        To turn multiple nodes into skip nodes, see also `skip_remaining`.
        """

        node = self.node_data(node_id)

        if node["expanded"]:
            return False

        pn = self.node_percolated_petri_net(node_id, compute=True)
        minimal_traps = trappist(network=pn, problem="min")
        minimal_traps = [(node["space"] | x) for x in minimal_traps]

        if len(minimal_traps) == 1 and minimal_traps[0] == node["space"]:
            # This node is a minimal trap space
            # and thus cannot be skipped.
            node["expanded"] = True
            return True

        for m_trap in minimal_traps:
            m_id = self._ensure_node(node_id, m_trap)
            # Also expand the minimal trap space, since we know
            # it has no successors.
            m_data = self.node_data(m_id)
            m_data["expanded"] = True

        node["expanded"] = True
        node["skipped"] = True

        if self.config["debug"]:
            print(f"[{node_id}] Added {len(minimal_traps)} skip edges.")

        return True

    def skip_remaining(self) -> int:
        """
        Apply `skip_to_minimal` to every node that is not expanded.

        This is faster than calling the method individually if the number of
        nodes is high since we can cache the minimal trap spaces.

        Returns the number of created skip nodes.
        """

        pn = self.node_percolated_network(self.root(), compute=True)
        root_space = self.node_data(self.root())["space"]
        minimal_traps = trappist(network=pn, problem="min")
        minimal_traps = [root_space | x for x in minimal_traps]

        if self.config["debug"]:
            print(f"Skipping remaining nodes. Found {len(minimal_traps)} trap spaces.")

        trap_with_id: list[tuple[int, BooleanSpace]] = []
        for m_trap in minimal_traps:
            m_id = self._ensure_node(None, m_trap)
            m_data = self.node_data(m_id)
            trap_with_id.append((m_id, m_trap))
            # We can expand the minimal trap spaces since
            # they don't have successors.
            m_data["expanded"] = True

        skipped_nodes = 0
        for node_id in self.node_ids():
            node = self.node_data(node_id)

            if node["expanded"]:
                continue

            skip_edges = 0
            for m_id, m_trap in trap_with_id:
                if is_subspace(m_trap, node["space"]):
                    self._ensure_edge(node_id, m_id, m_trap)
                    skip_edges += 1

            node["skipped"] = True
            node["expanded"] = True
            skipped_nodes += 1

            # At this point, all minimal traps must be expanded,
            # hence we should never skip one.
            assert not self.node_is_minimal(node_id)

        if self.config["debug"]:
            print(f"Skipped {skipped_nodes} nodes.")

        return skipped_nodes

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

        # If the node had any attractor data computed as unexpanded, these are
        # no longer valid and need to be erased.
        node["attractor_seeds"] = None
        node["attractor_candidates"] = None
        node["attractor_sets"] = None

        current_space = node["space"]

        if self.config["debug"]:
            print(
                f"[{node_id}] Expanding: {len(self.node_data(node_id)['space'])} fixed vars."
            )

        if len(current_space) == self.network.variable_count():
            # This node is a fixed-point. Trappist would just
            # return this fixed-point again. No need to continue.
            if self.config["debug"]:
                print(f"[{node_id}] Found fixed-point: {current_space}.")
            node["expanded"] = True
            return

        # We use the non-propagated Petri net for backwards-compatibility reasons here.
        # The SD created from the restricted Petri net is technically correct, but can
        # propagate some of the input values further and yields a smaller SD.
        source_nodes = []
        if node_id == self.root():
            source_nodes = extract_source_variables(self.petri_net)

        sub_spaces: list[BooleanSpace]

        # Only use the percolated PN if it is already known.
        pn = node["percolated_petri_net"]
        if pn is not None:
            # We have a pre-propagated PN for this sub-space, hence we can use
            # that to compute the trap spaces.
            partial_sub_spaces = trappist(
                pn,
                problem="max",
                optimize_source_variables=source_nodes,
                solution_limit=self.config["max_motifs_per_node"],
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
                solution_limit=self.config["max_motifs_per_node"],
            )

        # Release the Petri net once the sub_spaces are computed.
        # It might be needed later for attractor computation, but it
        # uses a lot of memory in large diagrams to keep all the nets
        # in memory.
        node["percolated_petri_net"] = None

        if len(sub_spaces) == self.config["max_motifs_per_node"]:
            raise RuntimeError(
                f"Exceeded the maximum amount of stable motifs per node ({self.config['max_motifs_per_node']}; see `SuccessionDiagramConfiguration.max_motifs_per_node`)."
            )

        # Sort the spaces based on a unique key in case trappist is not always
        # sorted deterministically.
        sub_spaces = sorted(
            sub_spaces, key=lambda space: space_unique_key(space, self.network)
        )

        if len(sub_spaces) == 0:
            if self.config["debug"]:
                print(f"[{node_id}] Found minimum trap space: {current_space}.")
            node["expanded"] = True
            return

        if self.config["debug"]:
            print(f"[{node_id}] Found sub-spaces: {len(sub_spaces)}")

        for sub_space in sub_spaces:
            self._ensure_node(node_id, sub_space)

        # If everything else worked out, we can mark the node as expanded.
        node["expanded"] = True

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
                percolated_network=None,
                percolated_petri_net=None,
                percolated_nfvs=None,
                attractor_candidates=None,
                attractor_seeds=None,
                attractor_sets=None,
                parent_node=parent_id,
                skipped=None,
            )
            self.node_indices[key] = child_id
        else:
            child_id = self.node_indices[key]

        assert child_id is not None

        if parent_id is not None:
            self._ensure_edge(parent_id, child_id, stable_motif)

        return child_id

    def _ensure_edge(self, parent_id: int, child_id: int, stable_motif: BooleanSpace):
        # TODO: It seems that there are some networks where the same child
        # can be reached through multiple stable motifs. Not sure how to
        # approach these... but this is probably good enough for now.
        if not self.dag.has_edge(parent_id, child_id):  # type: ignore
            self.dag.add_edge(parent_id, child_id, motif=stable_motif)  # type: ignore
        self._update_node_depth(child_id, parent_id)
