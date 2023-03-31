import os
import sys
import re

RE_SD_SIZE = re.compile("\\s*Succession diagram size: (\\d+)\\s*")
RE_STUB_NODES = re.compile("\\s*Stub nodes: (\\d+)\\s*")
RE_TRAPS = re.compile("\\s*Minimal traps: (\\d+)\\s*")
RE_ATTRACTORS = re.compile("\\s*Attractors: (\\d+)\\s*")
RE_TIME = re.compile("\\s*real\\s*(\\d+\\.?\\d*)\\s*")

result_dir = sys.argv[1]

for file in sorted(os.listdir(result_dir)):
    if not file.endswith("_out.txt"):
        continue

    # name, min. traps, attractors, SD size, stubs, time
    row = [file,0,0,0,0,0]

    with open(f"{result_dir}/{file}") as contents:
        for line in contents:
            if RE_TRAPS.match(line):
                row[1] = int(RE_TRAPS.match(line).group(1))
            elif RE_ATTRACTORS.match(line):
                row[2] = int(RE_ATTRACTORS.match(line).group(1))            
            elif RE_SD_SIZE.match(line):
                row[3] = int(RE_SD_SIZE.match(line).group(1))
            elif RE_STUB_NODES.match(line):
                row[4] = int(RE_STUB_NODES.match(line).group(1))
            elif RE_TIME.match(line):
                row[5] = round(float(RE_TIME.match(line).group(1)), 2)
                # The apostrophe prevents Google sheets from reading this as a date :(
                row[5] = f"'{row[5]}"

    print(", ".join([str(x) for x in row]))