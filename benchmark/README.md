# Benchmarking

This folder contains several scripts useful for benchmarking. 
All `bench_*` scripts solve one "benchmark problem" for a single given network file.
The last line of output of such computation is always a comma-separated list of values with relevant results.
You can then use `run_bench.py` to execute such benchmark script for all networks in a folder.
This will collect the full output of each script, as well as create an agggregate result table.
Each row of such table contains the network name, runtime of the whole script, and its output.

## Install `balm`

The benchmark scripts assume `balm` is installed as a package. 
To do this, go to the repository root and run `pip install .` 
(this step reqires at least `pip` version `23`).

## Example

For example, to test how long it takes to generate the succession diagram for minimal trap spaces,
you can use the following:

```
# 1h - timeout
# ../bbm-bnet-inputs-true - folder with networks
# bench_sd_expand_min.py - benchmark script
python3 run_bench.py 1h ../bbm-bnet-inputs-true bench_sd_expand_min.py
```

You can optionally add `-p X` or `-i` flag at the end of the argument list. 
The `-p` flag causes the benchmark to run in parallel in `X` processes.
The `-i` flag causes the benchmark to run in "interactive" mode where the user can skip
individual benchmarks.