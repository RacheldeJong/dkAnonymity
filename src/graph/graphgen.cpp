#include "graphgen.h"

sparsegraph get_empty_graph(const int n){
   int m;
   DYNALLSTAT(sparsegraph, g1, g1_sz);
   m = SETWORDSNEEDED(n);
   DYNALLOC2(sparsegraph, g1, g1_sz, n, m, "malloc");
   EMPTYGRAPH(g1, m, n);
   g1->nv = n;
   return *g1;
}

sparsegraph random_graph(const int n, const int p1, const int p2){
   DYNALLSTAT(sparsegraph, g1, g1_sz);
   *g1 = get_empty_graph(n);
   rangraph2_sg(g1, FALSE, 1, 1000, n);
   return *g1;
}

sparsegraph read_graph_from_file(const char* file_name, const int n){
   FILE* file;
   SG_DECL(sg1);
   
   if((file = fopen(file_name, "r")) == NULL){ 
      printf("Error: File %s could not be opened\n", file_name);
      return sg1;
   }
   readgraph_sg(file, &sg1, 1, 0, 1, n);
   fclose(file);

   return sg1;
}

sparsegraph get_ingoing_graph(const sparsegraph sg1){
   SG_DECL(sg2);
   size_t pos, to, deg, max_e_id;
   char c;
   size_t i, j, k;
   
   SG_ALLOC(sg2, sg1.nv, sg1.nde,"malloc");
   sg2.nv = sg1.nv;
   sg2.nde = sg1.nde;
   for(i = 0; i < sg1.nv; i++){
      sg2.d[i] = 0;
   }
   
   // Get number of incomming edges for each node
   max_e_id = 0;
   for(i = 0; i < sg1.nv; i++){
      pos = sg1.v[i];

      // Edge id
      for(j = pos; j < pos + sg1.d[i]; j++){
         to = sg1.e[j];
         sg2.d[to] += 1;
         max_e_id = j;
      }
   }

   // Set all edges to -1 (unused)
   for(i = 0; i < sg2.nde; i++)
      sg2.e[i] = -1;

   // Set correct start position for all nodes
   sg2.v[0] = 0;
   for(i = 1; i < sg2.nv; i++){
      sg2.v[i] = sg2.v[i-1] + sg2.d[i-1];
   }

   // Fill with edges
   for(i = 0; i < sg1.nv; i++){
      pos = sg1.v[i]; // List of nodes starts here
      for(j = pos; j < pos + sg1.d[i]; j++){ // Iterate over edges of node i
         to = sg1.e[j];
         
         k = 0;
         // Find empty spot in edge array
         while(sg2.e[sg2.v[to] + k] != -1) k++;
         if(k >= sg2.d[to]) printf("Error target %ld. k= %ld. to= %ld. edges = %d. pos=%ld\n", to, k, i, sg2.d[to], sg2.v[to]);
         // Add reversed edge: i->to becomes to->i
         sg2.e[sg2.v[to] + k] = i;
      }
   }
   return sg2;
}
