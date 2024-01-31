from biodivine_aeon import BooleanNetwork, AsynchronousGraph, ColoredVertexSet, VariableId, BddVariable, UpdateFunction
from typing import Generator
import sys
import random
import copy

## Requires at least biodivine_aeon 1.0.0a3

def erase_inputs(model: BooleanNetwork) -> BooleanNetwork:
	"""
	Substitute all constant inputs (i.e. true/false) and identity inputs (i.e. update function
	is identity) for unknown "free" functions.
	"""
	graph = AsynchronousGraph(model)
	inputs: list[VariableId] = []
	for var in model.variables():
		# Do the semantic check on the update function BDD.
		fn_bdd = graph.mk_update_function(var)
		identity_fn = UpdateFunction.mk_var(model, var)
		identity_bdd = graph.symbolic_context().mk_update_function(identity_fn)
		if fn_bdd.is_true() or fn_bdd.is_false():
			inputs.append(var)
		if fn_bdd == identity_bdd:
			inputs.append(var)

	# Erase all input functions.
	for input in inputs:
		model.set_update_function(input, None)
	
	# The operation probably made some regulations unused. Get rid of them.
	return model.infer_valid_graph()


def sample_random_models(
	model: BooleanNetwork, 
	count: int = 32, 
	seed: int = 0
) -> Generator[tuple[list[bool], BooleanNetwork], None, None]:
	# Here, we are sampling across all constants (`true`/`false`) and inputs 
	# (free function, identity function). If you want to sample over a different
	# set of "inputs", skip the `erase_inputs` function, but make sure the inputs 
	# you want to sample have a free (unknown) update function and no regulators.
	model = erase_inputs(model.infer_valid_graph())
	stg = AsynchronousGraph(model)
	all_colors = stg.mk_unit_colors()			
	if all_colors.is_singleton():
		# This model has no inputs.
		yield ([], model)
	else:		
		ctx = stg.symbolic_context()
		bdd_vars = ctx.bdd_variable_set()
		rng = random.Random(seed)

		# Prepare a mapping from variable names to their corresponding symbolic 
		# parameter variables. We know that all parameters are just constants,
		# hence their function table must have only one element.
		input_symbolic_var: dict[VariableId, BddVariable] = {}
		for var in model.implicit_parameters():
			table = ctx.get_function_table(var)
			assert len(table) == 1
			symbolic_var = table[0][1]
			input_symbolic_var[var] = symbolic_var

		print(f" >> Sampling... ")
		for s in range(count):					
			if all_colors.is_empty():
				# This model has fewer input valuations than the sample count.
				break

			print(f"{s+1}; ", end="")

			# Here, we pick a random valuation and subtract it from the set of colors.
			# Unfortunately, there doesn't seem to be a nicer way to do this at the moment
			# with a deterministic seed that is tied to the Python rng.

			valuation_seed = rng.randrange(0, 2**30)
			valuation_sample = all_colors.to_bdd().valuation_random(valuation_seed)
			assert valuation_sample is not None
			valuation_bdd = bdd_vars.mk_conjunctive_clause(valuation_sample)
			valuation_set = ColoredVertexSet(ctx, valuation_bdd)
			all_colors = all_colors.minus(valuation_set.colors())

			signature: list[bool] = []
			for (var, bdd_var) in input_symbolic_var.items():
				if valuation_sample[bdd_var]:
					signature.append(True)
					model.set_update_function(var, "true")
				else:
					signature.append(False)
					model.set_update_function(var, "false")

			# Make a copy of the model to avoid future modifications.
			yield (signature, copy.copy(model))
		print(" >> Done.")

if __name__ == '__main__':
    bn = BooleanNetwork.from_file(sys.argv[1])
    count: int = int(sys.argv[2]) if len(sys.argv) >= 3 else 32
    seed: int = int(sys.argv[3]) if len(sys.argv) >= 4 else 0
    for (signature, model) in sample_random_models(bn, count, seed):
    	print([int(x) for x in signature], model)