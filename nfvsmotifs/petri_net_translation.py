from __future__ import annotations
"""
    Implements the translation from a `BooleanNetwork` object into a Petri net that can be processed
    by Trappist. The Petri net is represented as a `DiGraph`, with nodes having either a `kind=place`
    or a `kind=transition` attribute.

    The variable names in the network have to be "sanitized" before translation. In particular, this 
    means they can't use any special characters beyond "_". In the resulting Petri net, we then use 
    a "b0_*" and "b1_*" prefix to distinguish the "zero" and "one" places created for each variable.
    This is also important for `clingo`, as we have to guarantee that in our logic program, symbols 
    start with a lowercase letter.

    Right now, we are using implicant enumeration through PyEDA BDDs for translating the update 
    functions into a series of PN transitions. Nevertheless, recent tests seem to indicate that
    AEON BDDs are actually quite a bit faster than PyEDA BDDs. In the future, we might want to 
    use them instead of PyEDA.
"""

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from biodivine_aeon import BooleanNetwork # type: ignore    
    from pyeda.boolalg.bdd import BDDNode # type: ignore

from networkx import DiGraph # type: ignore
from pyeda.boolalg.bdd import expr2bdd, bddvar
from nfvsmotifs.pyeda_utils import aeon_to_pyeda 
import re

def sanitize_network_names(network: BooleanNetwork, check_only: bool = False):
    """
        Verifies that all network variable names contain only alphanumeric characters and "_".
        If this is not the case, attempts to rename the variables to make them compliant.

        Note that AEON should already prune away most of the special characters when parsing models,
        but it still allows names like `Gene_{Subscript}` (i.e. curly brackets) which we have
        to sanitize here. 

        If `check_only=True` is specified, no renaming takes place and the function fails with
        an `Exception` instead.
    """
    extra_id = 0 # We use extra id to resolve possible name clashes.
    for var in network.variables():
        extra_id += 1
        name = network.get_variable_name(var)
        if not re.match("^[a-zA-Z0-9_]+$", name):
            if check_only:
                raise Exception(f"Found unsanitized variable: `{name}`.")
            # Replace all invalid characters with an underscore
            new_name = re.sub("[^a-zA-Z0-9_]", "_", name)
            try:
                network.set_variable_name(var, new_name)
            except:
                # Name clash happened. Try to resolve it using a unique suffix.
                # In the unlikely event that this too fails, fail the whole function.                
                new_new_name = f"{new_name}_id{extra_id}"
                try:
                    network.set_variable_name(var, new_new_name)
                except:
                    raise Exception(f"Cannot sanitize variable `{name}`. Both {new_name} and {new_new_name} are invalid.")

    # Technically, network is modified in place, but we might as well
    # return it for convenience.
    return network

def variable_to_place(variable: str, positive: bool) -> str:
    """
        Convert the name of a network variable to the name of the corresponding positive 
        or negative Petri net place.
    """
    if positive:
        return f"b1_{variable}"
    else:
        return f"b0_{variable}"

def place_to_variable(place: str) -> tuple[str, bool]:
    """
        Convert the name of a Petri net place to the name of the network variable, plus
        a Boolean indicating whether the original place was positive or negative.
    """
    if place.startswith("b1_"):
        return (place[3:], True)
    elif place.startswith("b0_"):
        return (place[3:], False)
    else:
        raise Exception(f"Invalid place name: `{place}`.")

def extract_variable_names(encoded_network: DiGraph) -> list[str]:
    """
    Extract the variable names from a Petri net encoded Boolean network.
    
    The variables are  sorted lexicographically, since the original BN ordering is not 
    preserved by the Petri net. However, BNs order variables lexicographically by default,
    so unless the Petri net was created from a custom BN (i.e. not from a model file),
    the ordering should be the same.
    """
    variables = []
    for node in encoded_network.nodes():
        if node.startswith("b0_"):
            variables.append(place_to_variable(node)[0])
    
    return sorted(variables)

def network_to_petrinet(network: BooleanNetwork) -> DiGraph:
    """
        Converts a `BooleanNetwork` to a `DiGraph` representing a Petri net encoding
        of the original network. For details about the encoding, see module description.

        Note that the given network needs to have "sanitized" names, otherwise the 
        method will fail (see `sanitize_network_names` in this module).
    """
    # Assert that all network names are already sanitized.
    sanitize_network_names(network, check_only=True)

    assert network.num_parameters() == 0, "Unsupported: Network contains explicit parameters."
    assert network.num_implicit_parameters() == 0, "Unsupported: Network contains implicit parameters."

    pn = DiGraph()

    # Create a positive and a negative place for each variable.
    places = {}
    for var in network.variables():
        name = network.get_variable_name(var)
        p_name = variable_to_place(name, positive=True)
        n_name = variable_to_place(name, positive=False)
        pn.add_node(p_name, kind="place")
        pn.add_node(n_name, kind="place")
        places[name] = (n_name, p_name)

    # Create PN transitions for implicants of every BN transition.
    t_id = 0
    for var in network.variables():
        var_name = network.get_variable_name(var)
        function = network.get_update_function(var)
        function = aeon_to_pyeda(function)
        
        vx = bddvar(var_name)
        fx = expr2bdd(function)

        positive_implicants = fx & ~vx
        negative_implicants = ~fx & vx

        # Add 0->1 edges.
        _create_transitions(pn, places, var_name, positive_implicants, go_up=True)
        _create_transitions(pn, places, var_name, negative_implicants, go_up=False)

    return pn
        

def _create_transitions(
    pn: DiGraph, 
    places: dict[str, tuple[str, str]], 
    var_name: str, 
    implicant_bdd: BDDNode, 
    go_up: bool
):        
    """
        Just a helper method that creates PN transitions from BDDs representing
        positive/negative BN transitions.
    """
    dir_str = "up" if go_up else "down"
    for t_id, implicant in enumerate(implicant_bdd.satisfy_all()):                
        t_name = f"tr_{var_name}_{dir_str}_{t_id + 1}"
        pn.add_node(t_name, kind="transition")
        # The transition moves a token either from "zero place" to the 
        # "one place", or vice versa.
        if go_up:
            pn.add_edge(places[var_name][0], t_name)
            pn.add_edge(t_name, places[var_name][1])
        else:
            pn.add_edge(places[var_name][1], t_name)
            pn.add_edge(t_name, places[var_name][0])
        for variable, value in implicant.items():  
            variable = str(variable) # Convert from BDD variable to name.
            if variable == var_name:
                continue
            # For the remaining variables, we simply check if the required
            # token is present in the corresponding place.
            pn.add_edge(places[variable][value], t_name)
            pn.add_edge(t_name, places[variable][value])
