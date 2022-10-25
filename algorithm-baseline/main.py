"""
A place for prototyping attractor detection algorithms.

Basic terminology:

 - A space is a dictionary where variable names are its keys and values are strings "0"/"1".
 If a variable is free, it is not present in the dictionary.

 - Percolation refers to the syntactic process of propagating variables fixed in a space using the update
 functions of the specific BN.

"""
from biodivine_aeon import *
from percolation import percolate
from pathlib import Path
from pyeda.boolalg.expr import expr
from trappist import compute_trap_spaces
from conversions import aeon_to_petri_net, space_to_aeon_set
from bnet import read_bnet

import sys

model_path = Path(sys.argv[1])
network = BooleanNetwork.from_bnet(model_path.read_text())

space = {}
space['v_External_Activator'] = "0"

print("Space", space)
print("percolates to", percolate(network, space))

petri_net = aeon_to_petri_net(network)
print(network)
print(petri_net)

stg = SymbolicAsyncGraph(network)

for space in compute_trap_spaces(petri_net, max_output=10, computation="max"):
	print(space)
	print(space_to_aeon_set(stg, space))
