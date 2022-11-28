from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from biodivine_aeon import BooleanNetwork # type: ignore

from networkx import DiGraph # type: ignore

def infer_signed_interaction_graph(network: BooleanNetwork) -> DiGraph:
    """
        Takes an arbitrary `BooleanNetwork` and extracts a signed interacion graph 
        based on the *actual* dependencies of the network's update functions.

        In particular, this eliminates things like unused variables and it will also
        infer correct monotonicity of inputs. Importantly, this is performed using
        BDDs, so it is not merely a syntactic transformation (but general BDD 
        limitations apply).

        In the resulting digraph, we use edge attribute `sign="+"` or `sign="-"`
        to differentiate between positive and negative edges. An iteraction that
        is not monotonic is represented as `sign="?"`.
    """
    # Since the model files can be incorrect, we have to infer the regulations
    # from the network functions directly. Note that this step uses BDDs, so in the
    # unlikely event that the update functions are too complex to represent as BDDs,
    # this will crash on an out-of-memory error, or just run for a very long time.
    network = network.infer_regulatory_graph()

    # Convert AEON `RegulatoryGraph` to a generic `DiGraph` that we can use later.
    rg = network.graph()
    ig = DiGraph()

    for var in rg.variables():
        ig.add_node(rg.get_variable_name(var))

    for reg in rg.regulations():
        if not reg['observable']:
            # In general, a fully specified network should only contain 
            # observable (i.e. essential) regulations.
            raise Exception("Unreachable: You are using this on partically specified networks, aren't you?")
        source = rg.get_variable_name(reg['source'])
        target = rg.get_variable_name(reg['target'])
        sign = "?"        
        if 'monotonicity' in reg:
            # Monotonic regulation---add only one edge.
            if reg['monotonicity'] == 'activation':
                sign = "+"
            elif reg['monotonicity'] == 'inhibition':
                sign = "-"
            else:
                raise Exception(f"Unreachable: unknown monotonicity {reg['monotonicity']}.")
        ig.add_edge(source, target, sign=sign)
    return ig