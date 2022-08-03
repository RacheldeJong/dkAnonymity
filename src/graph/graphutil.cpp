#include "graphutil.h"

std::unordered_set<int> get_neighborhood_nodes(const sparsegraph sg, const int v, const int d, int &edges){
   std::queue< int > q; // Queue with nodes to be visited
   std::unordered_set<int> vis; //Visited nodes (unordered_set for quick access when in set)
   
   int cur, nid;
   int depth = 0;
   int time_to_increase_depth = 1;
   edges = 0;
   vis.clear();
   q.push(v);
   vis.insert(v);

   // Iterate over nodes to be visited
   while(!q.empty()){
        
      cur = q.front();
      // Keep track of when to increase depth of search
      if((time_to_increase_depth -= 1) == 0){
         if(++depth > d) break;
         time_to_increase_depth = q.size();
      }
      q.pop();

      // Iterate over edges of current node
      for(int it = 0; it < sg.d[cur]; ++it){
         nid = sg.e[sg.v[cur] + it];
         edges++;
         // Do not add to queue again if node is already added
         if(vis.find(nid) != vis.end()) continue;
         q.push(nid);
         vis.insert(nid);
      }
   }
   return vis;
}

std::unordered_set<int> get_neighborhood_nodes_directed(const sparsegraph sgo, const sparsegraph sgi, const int v, const int d, int &edges){
   std::queue< int > q; // Queue with nodes to be visited
   std::unordered_set<int> vis; //Visited nodes (unordered_set for quick access when in set)
   
   int cur, nid;
   int depth = 0;
   int time_to_increase_depth = 1;
   edges = 0;
   vis.clear();
   q.push(v);
   vis.insert(v);

   // Iterate over nodes to be visited
   while(!q.empty()){
        
      cur = q.front();
      // Keep track of when to increase depth of search
      if((time_to_increase_depth -= 1) == 0){
         if(++depth > d) break;
         time_to_increase_depth = q.size();
      }
      q.pop();

      // Iterate over outgoing edges of current node
      for(int it = 0; it < sgo.d[cur]; ++it){
         nid = sgo.e[sgo.v[cur] + it];
         edges++;
         // Do not add to queue again if node is already visited
         if(vis.find(nid) != vis.end()) continue;
         q.push(nid);
         vis.insert(nid);
      }

      // Iterate over ingoing edges of current node
      for(int it = 0; it < sgi.d[cur]; ++it){
         nid = sgi.e[sgi.v[cur] + it];
         edges++;
         // Do not add to queue again if node is already visited
         if(vis.find(nid) != vis.end()) continue;
         q.push(nid);
         vis.insert(nid);
      }
   }

   return vis;
}

std::unordered_set<int> get_neighborhood_nodes_directed_distribution(const sparsegraph sgo, const sparsegraph sgi, const int v, const int d, int &edges){
   std::queue< int > q; // Queue with nodes to be visited
   std::unordered_set<int> vis; //Visited nodes (unordered_set for quick access when in set)
   
   int cur, nid;
   int depth = 0;
   int time_to_increase_depth = 1;
   edges = 0;
   vis.clear();
   q.push(v);
   vis.insert(v);

   // Iterate over nodes to be visited
   while(!q.empty()){
        
      cur = q.front();
      // Keep track of when to increase depth of search
      if((time_to_increase_depth -= 1) == 0){
         if(++depth > d) break;
         time_to_increase_depth = q.size();
      }
      q.pop();

      // Iterate over edges of current node
      for(int it = 0; it < sgo.d[cur]; ++it){
         nid = sgo.e[sgo.v[cur] + it];
         edges++;
         // Do not add to queue again if node is already visited
         if(vis.find(nid) != vis.end()) continue;
         q.push(nid);
         vis.insert(nid);
      }

      for(int it = 0; it < sgi.d[cur]; ++it){
         nid = sgi.e[sgi.v[cur] + it];
         edges++;
         // Do not add to queue again if node is already visited
         if(vis.find(nid) != vis.end()) continue;
         q.push(nid);
         vis.insert(nid);
      }
   }
   return vis;
}

void get_neighborhood(const sparsegraph sg1, sparsegraph &sg2,int &v, const int d){
   std::unordered_map<int, int> map; // Maps each node in sg to new node in sg2
   int deg, pos, pos1, to, edges;

   // Get nodes in k-neighborhood
   std::unordered_set<int> nodes = get_neighborhood_nodes(sg1, v, d, edges);
   int n = nodes.size();

   // Initialize sub-graph
   sg2.nv = nodes.size(); // Set number of nodes
   sg2.nde = 0;

   // Initialize map sg->sg2
   n = 0;
   for(auto it : nodes){
      map[it] = n;
      n++;
   }

   pos = 0;
   // For each node in the subgraph
   for(auto it : nodes){
      deg = 0;
      sg2.v[map[it]] = pos; // Set position of node in sg2
      pos1 = sg1.v[it]; 
      
      // Iterate over edges from it in original graph
      for(int i = 0; i < sg1.d[it]; i++){
         to = sg1.e[pos1 + i];

         // Check if edge is contained in subgraph
         if(nodes.find(to) == nodes.end()) continue;

         // Otherwise add edge to subgraph
         sg2.nde++;
         sg2.e[pos + deg] = map[to];
         deg++;
      }

      sg2.d[map[it]] = deg;
      pos += deg;
   }
   v = map[v]; // Set v to corresponding node in new graph
}

void get_neighborhood_directed(const sparsegraph sgo, const sparsegraph sgi, sparsegraph &sg2, int &v, const int d){
   std::unordered_map<int, int> map; // Maps each node in sg to new node in sg2
   int deg, pos, pos1, to, edges;

   // Get nodes in k-neighborhood
   std::unordered_set<int> nodes = get_neighborhood_nodes_directed(sgo, sgi, v, d, edges);
   int n = nodes.size();

   // Initialize sub-graph
   sg2.nv = nodes.size(); // Set number of nodes
   sg2.nde = 0;

   // Initialize map sg->sg2
   n = 0;
   for(auto it : nodes){
      map[it] = n;
      n++;
   }

   pos = 0;
   // For each node in the subgraph
   for(auto it : nodes){
      deg = 0;
      sg2.v[map[it]] = pos; // Set position of node in sg2
      pos1 = sgo.v[it];
      
      // Iterate over edges from it in original graph
      for(int i = 0; i < sgo.d[it]; i++){
         to = sgo.e[pos1 + i];

         // Check if edge is contained in subgraph
         if(nodes.find(to) == nodes.end()) continue;

         // Otherwise add edge to subgraph
         sg2.nde++;
         sg2.e[pos + deg] = map[to];
         deg++;
      }

      sg2.d[map[it]] = deg;
      pos += deg;
   }
   v = map[v]; // Set v to corresponding node in new graph
}

std::map<int, size_t> get_degree_distribution(const sparsegraph sg){
   std::map<int, size_t> degrees;
   int degree;

   for(size_t i = 0; i < sg.nv; i++){
      degree = sg.d[i];
      
      if(degrees.find(degree) == degrees.end()) degrees.insert({degree, 1});
      else degrees[degree]++;
   }
   return degrees;
}

std::map<int, size_t> get_degree_distribution_directed(const sparsegraph sg){
   std::map<int, size_t> outdegs;
   int outdeg;

   for(size_t i = 0; i < sg.nv; i++){
      outdeg = sg.d[i];
      
      if(outdegs.find(outdeg) == outdegs.end()) outdegs.insert({outdeg, 1});
      else outdegs[outdeg]++;
   }
   return outdegs;
}

long get_hash(sparsegraph sg){
   return hashgraph_sg(&sg, HASHKEY);
}

void print_graph(const sparsegraph sg){
   int n = sg.nv;
   int pos;

   for(int i = 0; i < n; i++){
      printf("%d: ", i);
      pos = sg.v[i];
      for(int j = 0; j < sg.d[i]; j++ ){
         printf("%d, ", sg.e[pos + j]);
      }
      printf("\n");
   }
}