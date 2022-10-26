from biodivine_aeon import *

def remove_static_constraints(network):	
	"""
		A method that removes all information about regulation monotonicity and
		essentiality from the given `BooleanNetwork`.

		This is mostly done to allow handling of randomly generated or otherwise
		machine pre-processed files that can contain subtle logical redundancies
		that AEON would otherwise detect as warnings.
	"""
	rg = RegulatoryGraph([network.get_variable_name(var) for var in network.variables()])
	for reg in network.graph().regulations():
		rg.add_regulation({
			'source': network.get_variable_name(reg['source']),
			'target': network.get_variable_name(reg['target']),
			'observable': False,
		})

	bn = BooleanNetwork(rg)
	for var in network.variables():
		bn.set_update_function(
			network.get_variable_name(var),
			network.get_update_function(var)
		)

	return bn

def has_parameters(network):
	if network.num_parameters() > 0:
		return True
	for var in network.variables():
		if network.get_update_function(var) == None:
			return True
	return False