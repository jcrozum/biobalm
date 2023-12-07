# type: ignore
import sys
import time

import pystablemotifs
from biodivine_aeon import BooleanNetwork

from nfvsmotifs._sd_algorithms.expand_source_SCCs import expand_source_SCCs
from nfvsmotifs.control import succession_control
from nfvsmotifs.SuccessionDiagram import SuccessionDiagram

sys.setrecursionlimit(150_000)
path = sys.argv[1]
print(f"model: {path}")

# Import the rules in nfvsmotifs.
t0 = time.perf_counter()
bn = BooleanNetwork.from_file(path)
t_import_nfvsmotifs = time.perf_counter() - t0
print(f"{t_import_nfvsmotifs = :>.3f}")

# Import the rules in pystablemotifs.
t0 = time.perf_counter()
primes = pystablemotifs.format.import_primes(path)
t_import_pystablemotifs = time.perf_counter() - t0
print(f"{t_import_pystablemotifs = :>.3f}")

print()


# Build the succession diagram in nfvsmotifs
t0 = time.perf_counter()
sd_bfs = SuccessionDiagram(bn)
sd_bfs.expand_bfs()
t_succession_bfs_nfvsmotifs = time.perf_counter() - t0
print(f"{t_succession_bfs_nfvsmotifs = :>.3f}")

t0 = time.perf_counter()
sd_scc = SuccessionDiagram(bn)
expand_source_SCCs(sd_scc)
t_succession_scc_nfvsmotifs = time.perf_counter() - t0
print(f"{t_succession_scc_nfvsmotifs = :>.3f}")


# Build the succession diagram in pystablemotifs
t0 = time.perf_counter()
ar = pystablemotifs.AttractorRepertoire.from_primes(
    primes, MPBN_update=True
)  # MPBN_update just for faster testing
t_succession_pystablemotifs = time.perf_counter() - t0
print(f"{t_succession_pystablemotifs = :>.3f}")

print()

# select target space
id = sd_bfs.minimal_trap_spaces()[0]
target = sd_bfs.node_space(id)


# do the control in nfvsmotifs
t0 = time.perf_counter()
interventions_nfvsmotifs = succession_control(bn, target)
t_control_target_nfvsmotifs = time.perf_counter() - t0
print(f"{t_control_target_nfvsmotifs = :>.3f}")

# do the control in nfvsmotifs
t0 = time.perf_counter()
interventions_nfvsmotifs3 = succession_control(bn, target, succession_diagram=sd_scc)
t_control_scc_nfvsmotifs = time.perf_counter() - t0
print(f"{t_control_scc_nfvsmotifs = :>.3f}")

# do the control in nfvsmotifs
t0 = time.perf_counter()
interventions_nfvsmotifs2 = succession_control(bn, target, succession_diagram=sd_bfs)
t_control_bfs_nfvsmotifs = time.perf_counter() - t0
print(f"{t_control_bfs_nfvsmotifs = :>.3f}")

# now do it in pystablemotifs
t0 = time.perf_counter()
interventions_pystablemotifs = ar.reprogram_to_trap_spaces(target)
t_control_pystablemotifs = time.perf_counter() - t0
print(f"{t_control_pystablemotifs = :>.3f}")
