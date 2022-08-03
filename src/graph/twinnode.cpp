#include "twinnode.h"

int twin_nbs = 5; // Default value
int twinnode_count;

std::vector<int> find_twin_nodes(const sparsegraph sg, std::map<int, std::set<int>> &twin_node_map){
   std::map<std::set<int>, int>targets; // Maps set of neighbours to a node
   std::set<int> neighbour_set;
   std::vector<int> node_set;
   twinnode_count = 0;

   // For each node, check if it has a twin node / clone
   for(size_t i = 0; i < sg.nv; i++){
      neighbour_set.clear();

      // Skip node if it has more than max_nbs neighbours
      if(sg.d[i] > twin_nbs){
         node_set.push_back(i);
         continue;
      }

      // Get set of neighbours
      for(size_t j = 0; j < sg.d[i]; j++)
         neighbour_set.insert(sg.e[sg.v[i] +j]);

      // Check if its neighbours are in targets list
      auto it = targets.find(neighbour_set);
      
      // If yes, it has a twin node
      if(it != targets.end()){

         (twin_node_map[it->second]).insert(i);

         twinnode_count++;
         continue;
      }
      
      // Otherwise, not: update target list and insert into all set
      else{
         targets[neighbour_set] = i;
         node_set.push_back(i);
      }
      neighbour_set.clear();
   }
   
   return node_set;
}

std::vector< std::vector<int> >process_twin_nodes(std::vector< std::vector<int> >eq, const std::map<int, std::set<int>>twin_node_map){
   std::vector< std::vector<int> >eq_new;
   std::vector<int> eq_class;
   std::set<int> nodes;

   // Iterate over classes
   for(auto it : eq){
      // For each class, add twin nodes
      for(auto it2 : it){
         eq_class.push_back(it2);

         // Check if node has a twin
         auto it3 = twin_node_map.find(it2);
         // Copy twins to new equivalence class
         if(it3 != twin_node_map.end()){
            std::copy(it3->second.begin(), it3->second.end(), std::back_inserter(eq_class));
         }
      }
      eq_new.push_back(eq_class);
      eq_class.clear();
   }
   return eq_new;
}