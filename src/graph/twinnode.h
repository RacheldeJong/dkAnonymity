// twinnode.h: 
//
// Contains code to find twin nodes and update equivalence classes
//
// Last edited: 3-8-2021
//

#ifndef TWINNODE
#define TWINNODE

#include "../util.h"
#include "../graph/graphutil.h"

extern int twin_nbs;
extern int twinnode_count;

// Find all twin nodes in sg and store them in twin_node_map
// e.g. if 1, 2, 3 are twins: (1, [2, 3]) is stored in twin_node_map
//  1 is the node that represents the twin nodes
//
// @sg1: Sparsegraph
// @twin_node_map: after execution, will contain the twin nodes of sg
//
std::vector<int> find_twin_nodes(const sparsegraph sg, std::map<int, std::set<int>> &twin_node_map);

// Updates eq by adding the twin nodes: the node representing a set of twins should be in one class of eq
// other twins can be added to this equivalence class
//
// @eq: the equivalence class without twin nodes
// @twin_node_map: the map created by twin_node_check that contains twin nodes
//
// @Returns: eq, where twin nodes are added in the correct class
//
std::vector< std::vector<int> >process_twin_nodes(std::vector< std::vector<int> >eq, const std::map<int, std::set<int>>twin_node_map);

#endif