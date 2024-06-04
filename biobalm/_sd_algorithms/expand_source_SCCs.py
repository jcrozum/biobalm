from __future__ import annotations

import itertools as it
from typing import TYPE_CHECKING, Callable, cast

from biobalm.space_utils import percolate_network
from biobalm.types import BooleanSpace
from biobalm.interaction_graph_utils import source_nodes

if TYPE_CHECKING:
    from biobalm.succession_diagram import SuccessionDiagram

    ExpanderFunctionType = Callable[[SuccessionDiagram], bool]


def expand_source_SCCs(
    sd: SuccessionDiagram,
    check_maa: bool,
    recursion: int = 0,
    expander: ExpanderFunctionType | None = None,
) -> bool:
    """
    1. percolate
    2. find source nodes, fix combinations and then percolate
    3. at every level, find all source SCCs, expand them
    4. when there are no more SCCs, expand in a usual way.

    TODO: it will be rare especially for empirical models,
    but there can be source SCCs even before fixing source node combinations.
    It could be useful to find them first, calculate scc_sd and store them somewhere,
    rather than calculating the same thing for every source node combination

    Parameters
    ----------
    check_maa - if True, check for motif avoidant attractors when
                expanding the sd with source SCCs.
                This grants significant speedup when searching for all attractors.
                if False, ignores motif avoidant attractors in source SCCs.

    """

    if expander is None:
        # By default, we recursively call the SCC expansion.
        def default_expander(sd: SuccessionDiagram):
            return expand_source_SCCs(sd, check_maa, recursion + 1)

        expander = default_expander

    if sd.config["debug"]:
        print(
            f"Start SD expansion using SCC decomposition on {sd.network}, recursion level {recursion}."
        )

    root = sd.root()

    # Usage of sets prevents node repetition in levels (this can happen if we independently
    # percolate to the same "downstream" node after fixing different source SCC values).
    current_level: set[int] = set([root])
    next_level: set[int] = set()

    # This already accounts for constant percolation.
    node_space = sd.node_data(root)["space"]

    # find source nodes
    perc_bn = percolate_network(sd.network, node_space)
    sources = source_nodes(perc_bn)

    if sd.config["debug"]:
        print(f" > Computed source/input variable(s): {sources}")

    # get source nodes combinations and expand root node
    if len(sources) != 0:
        # If there are too many source nodes, this can generate an absurdly large SD.
        # This would be a problem even without the SCC expansion, but we can just
        # stop the whole thing faster because we know how many nodes it generates now.
        if 2 ** len(sources) > sd.config["max_motifs_per_node"]:
            raise RuntimeError(
                f"Exceeded the maximum amount of stable motifs per node ({sd.config['max_motifs_per_node']}; see `SuccessionDiagramConfiguration.max_motifs_per_node`)."
            )
        else:
            if sd.config["debug"]:
                print(
                    f" > Expanding {len(sources)} source node into {2 ** len(sources)} SD nodes."
                )

        bin_values_iter = it.product(range(2), repeat=len(sources))
        for bin_values in bin_values_iter:
            valuation = cast(BooleanSpace, dict(zip(sources, bin_values)))
            sub_space = node_space | valuation

            next_level.add(sd._ensure_node(root, sub_space))  # type: ignore

        # This makes the root artificially "expanded". Also, there
        # can be no attractors here because we are just fixing the source nodes.
        sd.node_data(root)["expanded"] = True
        sd.node_data(root)["attractor_seeds"] = []
        sd.node_data(root)["attractor_sets"] = []
        current_level = next_level
        next_level = set()

    while len(current_level) > 0:
        if sd.config["debug"]:
            print(
                f" > Start SCC expansion of a BFS level with {len(current_level)} node(s)."
            )

        # For each node in the current level, we expand all source SCCs and put the
        # results into a new level.
        for node_id in sorted(current_level):
            source_scc_diagrams = list(sd.source_scc_subdiagrams(node_id))
            if sd.config["debug"]:
                print(
                    f" > [{node_id}] Found {len(source_scc_diagrams)} sub-diagrams while expanding node."
                )

            # If there are no source SCCs, this node is a fixed-point and we can save it
            # into the last level (every other node has at least one source SCC).
            if len(source_scc_diagrams) == 0:
                if sd.config["debug"]:
                    print(f"[{node_id}] > No source SCCs found. Node is a fixed-point.")
                assert len(sd.node_successors(node_id, compute=True)) == 0
                continue

            # If there is only one source SCC, we can do a normal expansion to get to the next level,
            # because the attach sub-diagram process would just expand it a copy it into node_id.
            # Furthermore, if we want to recursively use SCC expansion in SCC expansion, we *have to*
            # do this, otherwise we would just recursively call ourselves forver on a minimal trap spaces,
            # since the recursion wouldn't know where to stop.
            if len(source_scc_diagrams) == 1:
                if sd.config["debug"]:
                    print(f"[{node_id}] > Singe source SCCs found. Expanding normally.")
                next_level = next_level | set(sd.node_successors(node_id, compute=True))
                continue

            attach_at_list: list[int] = [node_id]
            for scc_diagram in source_scc_diagrams:
                fully_expanded = expander(scc_diagram)
                if not fully_expanded:
                    # Something bad happened in the expander function and we can't continue.
                    return False

                if sd.config["debug"]:
                    print(
                        f"[{node_id}] > Source SCC diagram expanded to {len(scc_diagram)} nodes."
                    )

                # At this point, diagram is fully expanded and we can attach its
                # nodes as the successors of `node_id`.
                next_attach_at_list: list[int] = []
                for attach_at in attach_at_list:
                    next_attach_at_list += attach_scc_subdiagram(
                        sd, scc_diagram, attach_at, check_maa
                    )
                attach_at_list = next_attach_at_list

            if attach_at_list == [node_id]:
                # Nothing was attached at this level. This means that all source SCCs
                # have a trivial succession diagram where we can't further "fix" anything else.
                # Usually, this means that we are in a minimal trap space. However, there can also
                # be another SCC connected to the source SCCs that can still have some
                # non-trivial trap space structure.
                #
                # Ideally, we should be able to "ignore" the current source SCCs and move on to the
                # following level of SCCs (because those source SCC stay "unstable" forever), but
                # that is currently not possible, so we just expand normally.
                if sd.config["debug"]:
                    print(
                        f"[{node_id}] > No further nodes attached. Expanding normally."
                    )
                next_level = next_level | set(sd.node_successors(node_id, compute=True))
            else:
                # Otherwise, move everything into the next layer.
                next_level = next_level | set(attach_at_list)

        current_level = next_level
        next_level = set()

    if sd.config["debug"]:
        print(
            f" > SCC expansion terminated with {len(sd)} node(s) on recursion level {recursion}."
        )

    return True


def attach_scc_subdiagram(
    sd: SuccessionDiagram, scc_sd: SuccessionDiagram, attach_at: int, check_maa: bool
) -> list[int]:
    """
    Attach the `scc_sd` to the `attach_at` node of `sd`.

    This means that every node of `scc_sd` will be extended to the full context of `sd`
    and attached as if `attach_at` was the root of `scc_sd`.

    If `check_maa` is set to `True`, the method will also check `scc_sd` for motif-avoidant
    attractors and if MAAs are disproven, the corresponding nodes of `sd` will be also marked
    as not having any MAAs.

    Parameters
    ----------
    sd - The main succession diagram into which we are attaching the `scc_sd`.
    scc_sd - A component subdiagram which will be copied into `sd`.
    attach_at - A node id from `sd` that will be used as the virtual "root" of `scc_sd` while attaching.
    check_maa - A flag that indicates that MAAs absence information should be propagated from the `scc_sd`
                to the main `sd`.

    Returns
    -------
    The IDs of `sd` nodes that correspond to the minimal trap spaces of `scc_sd`.
    """

    if len(scc_sd) == 1:
        # This `scc_sd` has a single minimal SCC, which thus corresponds to the `attach_at` point.
        return [attach_at]

    # Maps node IDs from the `scc_sd` to the extended and copied nodes in `sd`.
    node_id_map: dict[int, int] = {scc_sd.root(): attach_at}

    attach_at_space = sd.node_data(attach_at)["space"]

    # First, copy every node from `scc_sd`` into main `sd`.
    min_traps: list[int] = []
    for scc_node_id in scc_sd.node_ids():
        if scc_node_id == scc_sd.root():
            continue  # Root is implicitly copied.

        scc_node_space = scc_sd.node_data(scc_node_id)["space"]
        extended_node_space = scc_node_space | attach_at_space

        main_node_id = sd._ensure_node(parent_id=None, stable_motif=extended_node_space)  # type: ignore
        node_id_map[scc_node_id] = main_node_id

        if scc_sd.node_is_minimal(scc_node_id):
            min_traps.append(main_node_id)
        else:
            # This node can be marked as expanded, because we know its successors.
            # We just need to add them in the for loop below.
            sd.node_data(main_node_id)["expanded"] = True

        if check_maa:
            if len(scc_sd.node_attractor_candidates(scc_node_id, compute=True)) == 0:
                sd.node_data(main_node_id)["attractor_seeds"] = []
                sd.node_data(main_node_id)["attractor_sets"] = []

    assert len(node_id_map) == len(scc_sd)

    # Then copy all the edges.
    for scc_node_id in scc_sd.node_ids():
        for scc_node_succ in scc_sd.node_successors(scc_node_id):
            inner_stable_motif = scc_sd.edge_stable_motif(scc_node_id, scc_node_succ)

            main_node_id = node_id_map[scc_node_id]
            main_succ_id = node_id_map[scc_node_succ]

            # This should not happen, because the source SCCs are independent.
            # But I don't have a formal proof.
            assert main_node_id != main_succ_id

            sd._ensure_edge(main_node_id, main_succ_id, inner_stable_motif)  # type: ignore

    # This makes the `attach_at` node expanded. We will not be adding new nodes to it later.
    sd.node_data(attach_at)["expanded"] = True
    # Finally, if we are checking for MAAs, we can do that for the root too:
    if check_maa:
        if len(scc_sd.node_attractor_candidates(scc_sd.root(), compute=True)) == 0:
            sd.node_data(attach_at)["attractor_seeds"] = []
            sd.node_data(attach_at)["attractor_sets"] = []

    return min_traps
