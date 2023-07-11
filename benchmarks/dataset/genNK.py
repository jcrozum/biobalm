#!python3

from sys import argv, stdin, stdout
from typing import Generator, IO, List, Set
from subprocess import CalledProcessError
import subprocess


network_size = [100, 150, 200, 250, 300, 350, 400]
k = 2
n_instance = 20

# .r file
r_file = open("N-K.r", "w")
r_file.write("library(BoolNet)\n")

for n in network_size:
  for i in range(1, n_instance + 1):
    bn_file = str(n) + "-" + str(k) + "-" + str(i)

    line_cmd = f"net <- generateRandomNKNetwork(n={n}, k={k}, topology=\"homogeneous\", linkage=\"lattice\", functionGeneration=\"uniform\");"
    save_cmd = f"saveNetwork(net, file=\"NK/{bn_file}.bnet\");"

    line_cmd += save_cmd + "\n"
    r_file.write(line_cmd)

r_file.close()

subprocess.call(['Rscript', 'N-K.r'], stdout=subprocess.DEVNULL)




