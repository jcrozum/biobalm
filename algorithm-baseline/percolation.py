import pyeda.boolalg.expr as expr

# Percolate the given space according to the rules of the given network.
#
# Note that only AND/OR/NOT are supported at the moment.
def percolate(network, space):
	result = {}
	for var in network.variables():
		var_name = network.get_variable_name(var)

		# Transfer currently fixed variables.
		if var_name in space:
			result[var_name] = space[var_name]

		expression = expr.expr(network.get_update_function(var).replace("!", "~"))
		value = partial_eval(expression, space)		
		if value != None:
			result[var_name] = value
	return result

# Performs partial evaluation of the given pyeda expression within the given 
# space.
#
# TODO: Not all operators are supported at the moment, only AND/OR/NOT work.
def partial_eval(function, space):	
	if type(function) == expr.Variable:		
		key = str(function)
		if key in space:
			return space[key]
		return None
	if type(function) == expr.Complement:
		# Not sure what exactly is the difference between complement and NotOp,
		# but complement seems to only apply to literals.
		key = str(function.inputs[0])
		if key in space:
			if space[key] == "1":
				return "0"
			if space[key] == "0":
				return "1"
		return None
	if type(function) == expr.NotOp:
		inner = partial_eval(function.x, space)
		if inner == "0":
			return "1"
		if inner == "1":
			return "0"
		return None
	if type(function) == expr.AndOp:
		left = partial_eval(function.xs[0], space)
		right = partial_eval(function.xs[1], space)
		if left == "0" or right == "0":
			return "0"
		if left == "1" and right == "1":
			return "1"
		return None
	if type(function) == expr.OrOp:
		left = partial_eval(function.xs[0], space)
		right = partial_eval(function.xs[1], space)
		if left == "1" or right == "1":
			return "1"
		if left == "0" and right == "0":
			return "0"
		return None
	raise Exception(f"Unsupported expression: {function}; {type(function)}.")