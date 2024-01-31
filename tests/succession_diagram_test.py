import unittest

from biodivine_aeon import BooleanNetwork, AsynchronousGraph, Attractors

import balm
import balm.SuccessionDiagram
from balm.SuccessionDiagram import SuccessionDiagram
from balm.types import BooleanSpace

# This just ensures that the debug outputs are a part of the test output.
balm.SuccessionDiagram.DEBUG = True


class SuccessionDiagramTest(unittest.TestCase):
    def test_succession_diagram_structure(self):
        bn = BooleanNetwork.from_bnet(
            """
            x1, x2
            x2, x1
            x3, !x3
            """
        )

        succession_diagram = SuccessionDiagram(bn)
        succession_diagram.expand_bfs()
        assert succession_diagram.dag.number_of_nodes() == 3
        assert succession_diagram.dag.number_of_edges() == 2  # type: ignore
        assert max(d["depth"] for _, d in succession_diagram.dag.nodes(data=True)) == 1  # type: ignore
        assert succession_diagram.depth() == 1
        assert (
            max(
                [
                    succession_diagram.node_depth(i)
                    for i in succession_diagram.node_ids()
                ]
            )
            == succession_diagram.depth()
        )
        assert sum(1 for _ in succession_diagram.node_ids()) == len(succession_diagram)
        assert sum(1 for _ in succession_diagram.expanded_ids()) == len(
            succession_diagram
        )
        assert len(succession_diagram.minimal_trap_spaces()) == 2
        assert succession_diagram.find_node({"x1": 1, "x2": 0}) is None

        succession_diagram_one = succession_diagram

        bn = BooleanNetwork.from_bnet(
            """
        a, b
        b, a
        c, a & c & d | b & !c | c & !d
        d, !a | d | c
        """
        )

        succession_diagram = SuccessionDiagram(bn)

        # Initially, nothing is expanded so this should cause an error.
        with self.assertRaises(KeyError):
            succession_diagram.node_successors(succession_diagram.root())

        # Also, attractors are initially unknown too.
        with self.assertRaises(KeyError):
            succession_diagram.node_attractor_seeds(succession_diagram.root())

        # Expand the root manually and check that iterators work correctly.
        succession_diagram.node_successors(succession_diagram.root(), compute=True)
        assert (
            sum(1 for _ in succession_diagram.stub_ids()) == len(succession_diagram) - 1
        )
        assert sum(1 for _ in succession_diagram.expanded_ids()) == 1

        # Then expand the whole thing.
        succession_diagram.expand_bfs()
        assert succession_diagram.dag.number_of_nodes() == 4
        assert succession_diagram.dag.number_of_edges() == 5  # type: ignore
        assert max(d["depth"] for _, d in succession_diagram.dag.nodes(data=True)) == 2  # type: ignore
        assert succession_diagram.depth() == 2
        assert (
            max(
                [
                    succession_diagram.node_depth(i)
                    for i in succession_diagram.node_ids()
                ]
            )
            == succession_diagram.depth()
        )
        assert sum(1 for _ in succession_diagram.node_ids()) == len(succession_diagram)
        assert sum(1 for _ in succession_diagram.expanded_ids()) == len(
            succession_diagram
        )
        assert len(succession_diagram.minimal_trap_spaces()) == 2
        assert succession_diagram.find_node({"a": 1, "b": 0}) is None

        # The comparison functions should work even if the diagrams
        # are not based on the same network.
        assert not succession_diagram.is_subgraph(succession_diagram_one)
        assert not succession_diagram_one.is_subgraph(succession_diagram)
        assert not succession_diagram.is_isomorphic(succession_diagram_one)

        succession_diagram_partial = SuccessionDiagram(bn)
        succession_diagram_partial.node_successors(
            succession_diagram.root(), compute=True
        )

        assert succession_diagram_partial.is_subgraph(succession_diagram)
        assert not succession_diagram.is_subgraph(succession_diagram_partial)


def test_expansion_depth_limit_bfs():
    bn = BooleanNetwork.from_file("models/bbm-bnet-inputs-true/033.bnet")

    sd = SuccessionDiagram(bn)
    assert not sd.expand_bfs(bfs_level_limit=3)
    assert sd.expand_bfs(bfs_level_limit=10)
    assert len(sd) == 432


def test_expansion_depth_limit_dfs():
    bn = BooleanNetwork.from_file("models/bbm-bnet-inputs-true/033.bnet")

    sd = SuccessionDiagram(bn)
    assert not sd.expand_dfs(dfs_stack_limit=3)
    assert sd.expand_dfs(dfs_stack_limit=10)
    assert len(sd) == 432


def test_expansion_size_limit_bfs():
    bn = BooleanNetwork.from_file("models/bbm-bnet-inputs-true/033.bnet")

    sd = SuccessionDiagram(bn)
    assert not sd.expand_bfs(size_limit=200)
    assert sd.expand_bfs(size_limit=500)
    assert len(sd) == 432


def test_expansion_size_limit_dfs():
    bn = BooleanNetwork.from_file("models/bbm-bnet-inputs-true/033.bnet")

    sd = SuccessionDiagram(bn)
    assert not sd.expand_dfs(size_limit=200)
    assert sd.expand_dfs(size_limit=500)
    assert len(sd) == 432


# TODO: add tests for a wider variety of networks


def test_expansion_comparisons(network_file: str):
    # Compare the succession diagrams for various expansion methods.
    balm.SuccessionDiagram.DEBUG = True
    NODE_LIMIT = 100
    DEPTH_LIMIT = 10

    bn = BooleanNetwork.from_file(network_file)

    sd_bfs = SuccessionDiagram(bn)
    bfs_success = sd_bfs.expand_bfs(bfs_level_limit=DEPTH_LIMIT, size_limit=NODE_LIMIT)
    sd_dfs = SuccessionDiagram(bn)
    dfs_success = sd_dfs.expand_dfs(dfs_stack_limit=DEPTH_LIMIT, size_limit=NODE_LIMIT)

    if not (bfs_success and dfs_success):
        # succession_diagram is too large for this test.
        return

    sd_min = SuccessionDiagram(bn)
    assert sd_min.expand_minimal_spaces(size_limit=NODE_LIMIT)

    assert sd_bfs.is_isomorphic(sd_dfs)
    assert sd_min.is_subgraph(sd_bfs)
    assert sd_min.is_subgraph(sd_dfs)

    sd_attr = SuccessionDiagram(bn)
    assert sd_attr.expand_attractor_seeds(size_limit=NODE_LIMIT)

    assert sd_min.is_subgraph(sd_attr)
    assert sd_attr.is_subgraph(sd_bfs)

    # This will go through the minimal trap spaces of this network
    # and try to only expand towards this minimum trap space as a target.
    # This should always create a succession_diagram with exactly one minimal trap space,
    # as the rest
    for min_trap in sd_bfs.minimal_trap_spaces():
        space = sd_bfs.node_space(min_trap)

        sd_target = SuccessionDiagram(bn)
        assert sd_target.expand_to_target(space, size_limit=NODE_LIMIT)

        assert sd_target.is_subgraph(sd_bfs)
        assert len(sd_target.minimal_trap_spaces()) == 1


def test_attractor_detection(network_file: str):
    # TODO: Once attractor detection is faster, we should increase this limit.
    # Right now, checking attractors in larger succession diagrams would often
    # time out our CI.
    NODE_LIMIT = 100

    bn = BooleanNetwork.from_file(network_file)
    stg = AsynchronousGraph(bn.infer_valid_graph())

    # Compute the succession diagram.
    sd = SuccessionDiagram(bn)
    fully_expanded = sd.expand_bfs(bfs_level_limit=1000, size_limit=NODE_LIMIT)

    # succession_diagram must be fully expanded, otherwise we may miss some results.
    # If succession_diagram is not fully expanded, we just skip this network.
    if not fully_expanded:
        return

    # Compute attractors in diagram nodes.
    # TODO: There will probably be a method that does this in one "go".
    nfvs_attractors: list[BooleanSpace] = []
    for i in sd.node_ids():
        attr = sd.node_attractor_seeds(i, compute=True)
        for a in attr:
            # Just a simple sanity check.
            assert len(a) == bn.variable_count()
        if len(attr) > 0:
            nfvs_attractors += attr

    # Compute symbolic attractors using AEON.
    symbolic_attractors = Attractors.attractors(stg, stg.mk_unit_colored_vertices())

    # Check that every "seed" returned by SuccessionDiagram appears in
    # some symbolic attractor, and that every symbolic attractor contains
    # at most one such "seed" state.
    for seed in nfvs_attractors:
        symbolic_seed = stg.mk_subspace(seed)
        found = None

        # The "seed" state must have a symbolic attractor (and that
        # attractor mustn't have been removed yet).
        for i in range(len(symbolic_attractors)):
            if symbolic_seed.is_subset(symbolic_attractors[i]):
                found = i
        assert found is not None

        symbolic_attractors.pop(found)

    print("Attractors:", len(nfvs_attractors))

    # All symbolic attractors must be covered by some seed at this point.
    assert len(symbolic_attractors) == 0


def test_attractor_expansion(network_file: str):
    # This test is similar to the "test attractor detection" function above, but
    # it will perform only a partial expansion of the succession diagram, which
    # is hopefully faster.

    # TODO: Once attractor detection is faster, we should increase this limit.
    # Right now, checking attractors in larger succession diagrams would often
    # time out our CI.
    NODE_LIMIT = 100

    bn = BooleanNetwork.from_file(network_file)
    stg = AsynchronousGraph(bn.infer_valid_graph())

    # Compute the succession diagram.
    sd = SuccessionDiagram(bn)
    fully_expanded = sd.expand_attractor_seeds(size_limit=NODE_LIMIT)

    # succession_diagram must be fully expanded, otherwise we may miss some results.
    # If succession_diagram is not fully expanded, we just skip this network.
    if not fully_expanded:
        return

    # Compute attractors in diagram nodes.
    # TODO: There will probably be a method that does this in one "go".
    nfvs_attractors: list[BooleanSpace] = []
    # This is an important change compared to the original test: Here, we only
    # care about expanded nodes, everything else is ignored.
    for i in sd.expanded_ids():
        attr = sd.node_attractor_seeds(i, compute=True)
        for a in attr:
            # Just a simple sanity check.
            assert len(a) == bn.variable_count()
        if len(attr) > 0:
            nfvs_attractors += attr

    # Compute symbolic attractors using AEON.
    symbolic_attractors = Attractors.attractors(stg, stg.mk_unit_colored_vertices())

    # Check that every "seed" returned by SuccessionDiagram appears in
    # some symbolic attractor, and that every symbolic attractor contains
    # at most one such "seed" state.
    for seed in nfvs_attractors:
        symbolic_seed = stg.mk_subspace(seed)
        found = None

        # The "seed" state must have a symbolic attractor (and that
        # attractor mustn't have been removed yet).
        for i in range(len(symbolic_attractors)):
            if symbolic_seed.is_subset(symbolic_attractors[i]):
                found = i
        assert found is not None

        symbolic_attractors.pop(found)

    print("Attractors:", len(nfvs_attractors))

    # All symbolic attractors must be covered by some seed at this point.
    assert len(symbolic_attractors) == 0
