## Benchmarking

There is a "benchmark runner" written in Python in `benchmark/run_bench.py`. Currently, it works on linux and macOS (as long as you have `coreutils` installed). 
As arguments, you must give it (in this order):

 - Timeout string, as accepted by the `timeout` UNIX utility (e.g. `15s`, `1h`, etc.). Each benchmark will run with the given timeout.
 - A path to a directory with benchmark models. The script will automatically test all `.bnet`, `.sbml` and `.aeon` files in this directory.
 - A path to a binary which should be executed with each model in the benchmark directory (model is given as the first argument). If the path points to a Python script, benchmark runner will automatically start it using Python. Otherwise, it treats it as a native binary.
 - \[Optional] `-p NUM` to run the benchmarks in parallel on `NUM` cores. Note that running more benchmarks in parallel will typically degrade performance even if sufficient dedicated CPU cores are available.
 - \[Optional] `-i` to run in interactive mode. In interactive mode, the script will ask before starting the next benchmark, and each benchmark can be skipped. Not supported in parallel mode.

Benchmark runner will create a directory `_run_{BENCHMARK_DIRECTORY}_{BENCHMARK_BINARY}_{TIMESTAMP}`. In this directory, the output of each benchmark is stored, together with `.csv` files containing basic summary of runtimes.

**Limitations:**
 - The script *tries* to properly propagate and monitor for errors, but for some reason not all errors are always recognized correctly (e.g. out-of-memory tends to just fail silently). The only error that seems to propagate reliably is the timeout. It is recommended to check the output files for some easily verifiable signs of success (e.g. put a `Result: ...` line as the last output and then run a regex to verify that this is indeed part of the file).
 - For reasons that are unknown to me, Python only flushes its standard output in very long intervals. In case of a timeout, any results that are not written out are lost. In practice, this means that a benchmark which terminates due to a timeout can have non-trivial parts of its output missing. Flushing the output stream in your benchmark script helps, but it still does not seem to resolve this issue universally (I've only observed this behaviour when testing Python programs... other binaries seem to be fine as long as they flush regularly).

#### Examples

```
python run_bench.py 1h bbm-bnet-inputs-true bench_sd_construction.py -p 4
```

Run `bench_sd_construction.py` for all models in `bbm-bnet-inputs-true` with a timeout of one hour, running at most 4 processes in parallel.

```
python run_bench.py 45s PyBoolNet-repo bench_attractor_search.py
```

Run `bench_attractor_search.py` for all models in `PyBoolNet-repo` with a timeout of 45 seconds, a single model at a time.

