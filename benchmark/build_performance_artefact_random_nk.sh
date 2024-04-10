#!/bin/bash

# The purpose of this script is to compute an "artefact" that can be used 
# to compare the performance of `biobalm` between versions consistently.

# Note that even though we try to limit the evaluation as much as possible,
# the whole process will still typically require several hours.

# Also note that on laptops (and improperly built desktops :D), the
# temperature of the environment can have a measurable impact on performance
# and we are not doing any special variance analysis... always measure
# in a reasonably consistent conditions.

# The script will setup its own virtual environment in which Python
# dependencies will be installed. However, it still assumes you have pint
# and mole installed correctly.

# This script does not perform any correctness checking, but it could be 
# a good idea to validate that the results that were computed successfully
# actually match between versions.

if [[ $(ls _run_*) ]]; then
    echo "There are existing result runs already present. Please remove them first."
    exit 2
fi

if git diff-index --quiet HEAD --; then
    echo "No uncommitted changes. Installing biobalm..."
else
    echo "There are uncommitted changes. Please commit or stash them."
    exit 2
fi

rm -rf ./venv
python3 -m venv ./venv

# Install biobalm and dependencies.
./venv/bin/pip3 install ..

# Save the name of the machine and the current git commit for 
# future reference.
git rev-parse HEAD > _run_git_rev.txt
hostname > _run_hostname.txt

# Benchmark "basic" unconstrained expansion.
./venv/bin/python3 run_bench.py 10m ../models/random_nk bench_sd_expand_bfs.py
# Benchmark "basic" attractor enumeration.
./venv/bin/python3 run_bench.py 10m ../models/random_nk bench_sd_attractors_full.py
# Benchmark improved "SCC" attractor enumeration and expansion.
./venv/bin/python3 run_bench.py 10m ../models/random_nk bench_sd_attractors_scc.py
# Benchmark simple exhaustive control.
./venv/bin/python3 run_bench.py 10m ../models/random_nk bench_sd_control.py

zip -r perf_random_nk_`hostname`.zip _run_*