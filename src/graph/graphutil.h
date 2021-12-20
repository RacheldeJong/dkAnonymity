// graphutil.h
// 
// Contains utility functions required for computing the d-k-anonymity distribution
// of a graph
//
// Author: Rachel de Jong
// Date: 25-8-2020
//

#include "../util.h"

// Given a graph sg1 (sgo and sgi for directed), vertex v and distance d, 
// this function generates and returns the d-neighbourhood of v
// 
// @sg: sparse graph in which n-neighborhood should be found (undirected)
// @sgo: graph with outgoing edges only (directed)
// @sgi: graph with ingoing edges only (directed)
// @sg2: graph where d-neighbourhood is stored
// @v: node for which the d-neighborhood should be found
//     NOTE: this value changes based on the label it gets in the subgraph
//           nodes in the resulting graph (sg2) will have new labels
// @d: distance used for neighbourhood (d-neighbourhood)
//
void get_neighborhood(const sparsegraph sg, sparsegraph &sg2, int &v, const int d);
void get_neighborhood_directed(const sparsegraph sgo, const sparsegraph sgi, sparsegraph &sg2, int &v, const int d);

// Given a graph sg, this function finds all nodes that are in the d-neighborhood
// of node v. The number of edges it contains is also stored in variable edges
//
// @sg: sparse graph in which d-neighborhood should be found
// @v: node for which the n-neighborhood should be found
//     NOTE: this value changes based on the label it gets in the subgraph
// @k: k-neighborhood, maximum distance for nodes in subgraph from node vs
// @edges: Number of edges obtained in the subgraph. 
//         NOTE: This value is not used, but only changed during execution so the 
//         number of edges is known afterwards (Useful for allocating space for subgraph)
//
std::unordered_set<int> get_neighborhood_nodes(const sparsegraph sg, const int v, const int d, int &edges);

// Get the degree distribution (out-degree distribution for directed graphs)
// and store it in degrees (outdegs)
//
// @sg: sparse graph for which the degree distribution is computed
// @degrees / outdegs: data structure in which the degree distribution is stored
//
void get_degree_distribution(const sparsegraph sg, std::map<int, size_t>& degrees);
void get_degree_distribution_directed(const sparsegraph sg, std::map<int, size_t>& outdegs);

// Simple function to print sparse graph. Format:
// from: to, to,
// 0: 1, 2, 
// 1: 2, 3,
// ...
//
// @sg: sparse graph to print
// 
void print_graph(const sparsegraph sg);
