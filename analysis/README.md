## Example

```
# 1h - timeout
# ../bbm-bnet-inputs-true - folder with networks
# SD_analysis.py - benchmark script
python3 ../benchmark/run_bench.py 1h ../models/bbm-bnet-inputs-true SD_analysis.py
```

You can optionally add `-p X` or `-i` flag at the end of the argument list. 
The `-p` flag causes the benchmark to run in parallel in `X` processes.
The `-i` flag causes the benchmark to run in "interactive" mode where the user
can skip individual benchmarks.