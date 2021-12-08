# d-k-Anonymity
This repository contains code used in [1] to measure anonymity in complex networks.
Using this measure, vertices are partitioned into equivalence classes.
Vertices are equivalent if they have the same structural position in their d-neighbourhood; parameter d can be used to set the strictness based on different levels of how much a possible attacker knows.
This framework uses Nauty [2].

# The framework
* `src`: Directory containing source code
  * `anonymity.cpp`: Contains the main function and code
  * `util.h`: Includes the required libraries
  * Directory `graph`: Contains files related to graphs
      * `graphgen.cpp`: Code to generate graphs and read graphs from input
      * `graphutil.cpp`: Functions useful for working with graphs such as get_neighbourhood or get_degree_distribution
  * `Measure`: Contains files related to the graph measures
      * `dk-anonymity.cpp`: Code to partition vertices of a (possibly directed) graph based on d-k-anonymity
* `examples`: Directory containing example graphs

# Installation
Before using this code, the nauty framework should be downloaded from: https://pallini.di.uniroma1.it/
* Move this git-directory (dkAnonymity) in the downloaded nauty directory (nauty27r3)
* In the makefile in nauty27r3 add the following lines after line #500:


```
ANON = ./dkAnonymity/src
GRAPH = ${ANON}/graph
MEAS = ${ANON}/measure

anonymity :
			g++ -std=c++11 -o  ${ANON}/anonymity ${CFLAGS} ${ANON}/anonymity.cpp ${GRAPH}/graphutil.cpp ${GRAPH}/graphgen.cpp ${MEAS}/dk-anonymity.cpp traces.o nauty.a ${LDFLAGS}
```
`WARNING`: when cleaning the nauty directory, the makefile is regenerated and these lines are deleted.


For an older nauty version (nauty27r1) use: 
```
ANON = ./dkAnonymity/src
GRAPH = ${ANON}/graph
MEAS = ${ANON}/measure

anonymity : ${ANON}/anonymity.cpp
	        g++ -o ${ANON}/anonymity ${CFLAGS} ${ANON}/anonymity.cpp ${GRAPH}/graphutil.cpp ${GRAPH}/graphgen.cpp ${MEAS}/dk-anonymity.cpp traces.o nauty.a ${LDFLAGS}

```
`WARNING`: when cleaning the nauty directory, the makefile is regenerated and these lines are deleted.

# Compilation
```
make ./dkAnonymity/anonymity 
```

# How to run
To run the code, at least one input graph should be given. Command line arguments are optional.
Use the following command to run the code:
```
./dkAnonymity/src/anonymity ./examples/example2 [-arguments]
```
For example:
```
./dkAnonymity/src/anonymity ./dkAnonymity/examples/example2 -dir 1 -d 6
```

## Command line arguments
Various command line arguments can be used to adjust the settings of the algorithm.
* `-dir 0` : Set graph to directed or undirected
    * `0` - Undirected graphs (default)
    * `1` - Directed graphs
* `-d 5` : Set the distance (default 5)
* `-h 3` : Set the heuristic / algorithm choice
    * `0` - Naive (slowest)
    * `1` - Iterative
    * `2` - Iterative + degree distribution heuristic (fastest)
    * `3` - Iterative + number of vertices and edges as heuristic (fastest, default)
* `-eq 1`: Choose to print equivalence classes per iteration 
    * `1` - On (default)
    * `0` - Off
* `-s 1`: Set which statistics to print
    * `0` - Print final statistics only
    * `2` - Print statistics per iteration (default)
    * `3` - Print per class split (Preceded by `/`, debug mode)
    * `4` - Print all statistics, including progress updates (Preceded by `//`, debug mode)

# Input graphs

Input graph files should be using the .dre format. Some examples are included in directory 'examples'
Start with !n=number_vertices, then for each vertex 0 to n-1 vertex: edge1, target2, ...;. For the last vertex, end with `.` instead of `;`.

Example ` example2`:

```
!n=6
0: 3 4 5;
1: 3 4 5;
2: 3 4 5;
3: 0 1 2;
4: 0 1 2;
5: 0 1 2.
```
# Examples
Folder `examples` contains some small toy examples to test the d-k-anonymity measure. For each test case "x" there are 3 files:
* testx: The edgelist that can be read by Nauty
* testx.group: The grouping that d-k-anonymity will make if d = diameter(G). Vertices on the same line are in the same group
* testx.png: A visualization of the graph with each node coloured according to the group it belongs in. Generated with Networkx in `visualize.py`

Note that test1 is an empty graph: networkx does not show these nodes. Therefore `test1.png` is empty
Examples can be generated with `visualize.py`, code to generate the examples are given in `examples.sh`

# References
[1] B. D. McKay and A. Piperno, “Practical graph isomorphism”, Journal of Symbolic Computa-tion, vol. 60, no. 0, pp. 94–112, 2014. <br />
[2] A. Hagberg, D. A. Schult, and P. J. Swart, “Exploring network structure, dynamics, and functionusing networkx,” inProceedings of the 7th Python in Science Conference(G. Varoquaux, T. Vaught,and J. Millman, eds.), (Pasadena, CA USA), pp. 11 – 15, 2008. <br />
[3] G. Csardi and T. Nepusz, “The igraph software package for complex network research,”InterJournal,vol. Complex Systems, p. 1695, 2006. <br />
