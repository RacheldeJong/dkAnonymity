// dk-anonymity.h: 
// Contains functions to partition vertices of a graph with respect to 
// d-k-anonymity.
// 2 nodes are equivalent if:
// - Their d-neighborhoods are isomorph
// - There is an isomorphic function mapping the vertices onto each other
// Author: Rachel de Jong
// Data: 25-8-2021
//

#ifndef DKANONYMITY2
#define DKANONYMITY2

#include "../util.h"
#include "../graph/graphutil.h"

extern int heuristic_choice;
extern int print_eq_class;
extern int print_statistics;

// Returns true if sg1 and sg2 are the same wrt nodes v1 and v2 i.e.:
// - sg1 and sg2 are isomorphic
// - nodes v1 and v2 are isomorphic in sg1 sg2
// Otherwise returns false
//
// @sg1, sg2: Sparsegraphs
// @v1, v2: Nodes in sg1 and sg2 respectively
//
bool are_same_sg(sparsegraph *sg1, sparsegraph *sg2, const int v1, const int v2);

// Same as are_same_sg but the canonical labeling of the graphs are already computed and the required parts are passed to this function)
// Returns true if cg1 and cg2 are the same wrt nodes v1 and v2 i.e.:
// - cg1 and cg2 are isomorphic
// - nodes v1 and v2 are isomorphic in cg1 and cg2
// Otherwise returns false
// @cg1, cg2: Canonical labeling of graphs
// @v1, v2: Nodes in sg1 and sg2 respectively
// @v1_pos: Position of v1 in cg2
// @lab1, lab2, orbits: values obtained via sparsenauty function
//
bool are_same_sg_can(sparsegraph *cg1, sparsegraph *cg2, const int v1, const int v2, const int v1_pos, int *lab1, int *lab2, int *orbits);

// Given a graph compute the equivalence class for the k-neighborhood
// Each equivalence class will contain in the end a set of nodes such that:
// - The k-neighborhoods of all nodes are isomorphic
// - The the nodes themselves are isomorphic in their k-neighborhoods
//
// @sg: Sparse graph
// @n: Number of nodes in the graph
// @d: Neighborhood distance
//
std::vector< std::vector< int > > get_equivalence_classes(sparsegraph sg, const int d);
std::vector< std::vector< int > > get_equivalence_classes_directed(sparsegraph sgo, sparsegraph sgi, const int d);

// Given an equivalence class, this is split in multiple equivalence classes looking at the d-neighborhood
// (This function is used in are_same_sg)
//
// @sg: Sparse graph
// @eclass: An equivalence class
// @d: Neighbourhood distance
//
std::vector< std::vector< int > > split_equivalence_class(sparsegraph sg, std::vector <int> eclass, const int d);
std::vector< std::vector< int > > split_equivalence_class_directed(sparsegraph sgo, sparsegraph sgi, std::vector <int> eclass, const int d);

// Prints the given equivalence class
//
// @ eclasses: The equivalence classess
//
void print_equivalence_classes(const std::vector< std::vector < int > > eclasses);

// Prints the given equivalence class and outputs it to the given file
//
// @eclasses: The equivalence classes
// @file_name: File
//
void print_equivalence_classes_to_file(const std::vector < std::vector< int > > eclasses, char * file_name);

// Get k: size of the smallest equivalence class
//
// @elass: the equivalence classes
//
int get_k(const std::vector< std::vector< int > > eclasses);

// Print statistics about equivalence class
//
// @eclasses: Equivalence partition
// @id: Iteration number
void print_statistics_eq(const std::vector <std::vector <int > > eclasses, const int id);

#endif