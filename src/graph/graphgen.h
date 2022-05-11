// graphgen.h: 
//
// Contains code to generate graphs in Nauty:
// - Empty graph (no edges)
// - Random graph
// - From file)
// 
// Last edited: 11-5-2021
//

#include "../util.h"

#define ISDIGIT(c) ((c) >= '0' && (c) <= '9')

// Returns an empty nauty graph with n nodes. To this graph edges could be added.
//
// @n: Number of nodes the empty graph should contain
//
// @Returns: sparse graph with n nodes and 0 edges
//
sparsegraph get_empty_graph(const int n);

// Generate a random graph with n nodes and probability of p1 / p2 to contain each 
// possible edge
//
// @n: Number of nodes in graph
// @p1: Used in formula to compute probability that an edge exists
// @p2: Used in formula to compute probability that an edge exists
//
// @Returns: random sparse graph containing n nodes and a number of edges determinded by p1 and p2
//
sparsegraph random_graph(const int n, const int p1, const int p2);

// Reads and generates graph from file "file name" containing n nodes
// NOTE: the function does not check if n is correct. Edges from or to nodes
// that are not in range 0-(n-1) are not added to the graph
// 
// @file_name: String containing the path to a file containing the graph
// @n: Number of nodes the graph has
//
// @Returns: sparse graph representing graph stored in file_name
//
sparsegraph read_graph_from_file(const char* file_name, const int n);

// Get directed graph with reversed edges: incomming edges
// become outgoing edges
//
// @sg1: sparse graph
// 
// @Return: the resulting sparse graph where edges are reversed
// (required for directed variant of d-k-anonymity)
sparsegraph get_ingoing_graph(const sparsegraph sg1);
