

# Trappist: a tool for computing min./max. trap spaces and fixed points of Boolean networks 

## Requirements

+ Python 3
+ networkx
+ pyeda

## Usage

### Compute max. trap spaces

Run

    python3 trappist.py -c max -m 1000 -t 120 -s asp test.bnet

where 
+ `-c max` indicates computing max. trap spaces.
+ `-m 1000` indicates that the maximum number of solutions computed is `1000`.
+ `-t 120` indicates that the time limit is 120 seconds.
+ `asp` indicates that the problem is encoded as an ASP (it is possible to use MaxSAT or CSP).
+ `test.bnet` indicates the input file. Now, Trappist only supports the .bnet format.

Note that a Boolean network may have no max. trap spaces.

### Compute min. trap spaces

Run 

    python3 trappist.py -c min -m 1000 -t 120 -s asp test.bnet

where
+ `-c min` indicates computing min. trap spaces.

Note that a Boolean network always has at least one min. trap spaces.

### Compute fixed points

Run

​		python3 trappist.py -c fix -m 1000 -t 120 -s asp test.bnet

where
+ `-c fix` indicates computing fixed points.

Note that a Boolean network may have no fixed points.
