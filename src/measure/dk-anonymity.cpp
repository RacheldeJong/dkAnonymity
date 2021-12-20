#include "dk-anonymity.h"

int print_eq_class = 0;
int print_statistics = 1;
int heuristic_choice = 0;
int do_twin_node_check = 1;
int twin_node_count = 0;
int print_time_can_labelling = 0;
int twin_nbs = 3;
size_t iso_checks, nr_can1;

bool are_same_sg_can(sparsegraph *cg1, sparsegraph *cg2, const int v1, const int v2, const int v1_pos, int *lab2, int *orbits){
   int i, n, m, target;
   n = cg1->nv;
   m = SETWORDSNEEDED(n);

   iso_checks++;

   // Compare canonically labelled graphs
   // If graphs are isomorph, then after Traces function sg1 == sg2 holds
   if (aresame_sg(cg1,cg2)){
      
      // Find node in sg2 isomorph with v1:
      // Each node i in the graph is mapped to a new label: lab1[i](G1) -> i (G'1)
      // lab1[i] =iso lab2[i] -> if v1 = lab1[i] then lab2[i] is a node in G2 isomorphic to v1 in G1
      target = lab2[v1_pos];

      // Check if nodes are in the same orbit 
      // v1 =iso lab1[i] =iso lab2[i] and lab2[i] == target
      // -> v1 =iso target
      // if orbits[target] == orbits[v2] then target =iso v2 thus v1 =iso target =iso v2
      // -> v1 =iso v2
      if(orbits[v2] == orbits[target]){
         return true;
      }
   }
   return false;
}

std::vector< std::vector< int > > split_equivalence_class(sparsegraph sg, std::vector <int> eclass, const int d){
   std::vector< std::vector< int > >neweclass = {};
   std::vector< int > it_classes;
   std::vector<int> class_id(sg.nv, -1); // Maps node to class id
   SG_DECL(sub1); SG_DECL(sub2);
   SG_DECL(cansub1); SG_DECL(cansub2);
   int cur1, cur2, n, m, nodes, edges, i, v1_pos;
   int temp_val;
   size_t tel;
   bool added;
   clock_t t1;

   // Keys
   std::map< std::pair<int, int>, std::vector< int > > class_key; // Use nr. nodes + edges as map
   std::map< std::map<int, size_t>, std::vector< int > > degree_keys; // Use outdegree distribution as map
   std::map<int, size_t> degrees;

   // Declare required arrays
   DYNALLSTAT(int, lab, lab_sz);
   DYNALLSTAT(int, ptn, ptn_sz);
   DYNALLSTAT(int, orbits, orbits_sz);
   DYNALLSTAT(int, map, map_sz);
   static DEFAULTOPTIONS_SPARSEGRAPH(options);
   options.getcanon = TRUE;
   statsblk stats;
   
   // Return if class has size 1, can not split further
   if(d == 0 || eclass.size() == 1){
      neweclass.push_back(eclass);
      return neweclass;
   }
   n = sg.nv;
   m = sg.elen;

   // Allocate space for neighborhood graphs
   SG_ALLOC(sub1, n, m, "malloc");
   SG_ALLOC(sub2, n, m, "malloc");

   tel = 0;
   // Iterate over nodes in current eclass, place them in new subclass
   for(auto it : eclass){
      added = false;
      cur1 = it;

      //Update statistics:
      if(print_statistics >= 4 && (tel == 0 || tel % int((eclass.size() / 100) + 1) == 0)){
         printf("//node %ld / %ld\n", tel + 1, eclass.size());
         printf("//can nbh1 %ld\n", nr_can1);
         printf("//iso checks %ld\n", iso_checks);
         printf("//new classes %ld\n", neweclass.size());
         printf("//heuristic classes ");
         if(heuristic_choice == 2) printf("%ld\n", degree_keys.size());
         else printf("%ld\n", class_key.size());
         fflush(stdout);
      }
      tel += 1;

      // Get neighborhood of node to be placed
      get_neighborhood(sg, sub1, cur1, d);
      n = sub1.nv;
      m = sub1.nde;

      // Dynamic allocation for Traces
      DYNALLOC1(int, lab, lab_sz, n, "malloc");
      DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
      DYNALLOC1(int, orbits, orbits_sz,n, "malloc");
      DYNALLOC1(int, map, map_sz, n, "malloc");

      // Get keys using selected heuristic
      if(heuristic_choice == 1) // No heuristic
         it_classes = class_key[std::make_pair(0, 0)];
      else if(heuristic_choice == 2){ // Outdegree distribution
         get_degree_distribution(sub1, degrees);
         it_classes = degree_keys[degrees];
      }
      else // Use nr. nodes and edges of neighborhood (default)
         it_classes = class_key[std::make_pair(n, m)];

      // Try to place node in one of the found classes
      if(it_classes.size() != 0){

         // Get canonically labeled graph and position of it in lab1 array
         t1 = clock();
         sparsenauty(&sub1, lab, ptn, orbits, &options, &stats, &cansub1);
         if(print_time_can_labelling)
            printf("Can time first: %d, %f\n", sub1.nv, ((double)(clock() - (double)t1))/CLOCKS_PER_SEC);

         for (i = 0; i < n; i++){
            if(lab[i] == cur1) {
               v1_pos = i; 
               break;
            }
         }
         nr_can1++;

         for(int it2 = 0; it2 < it_classes.size(); it2++){
            
            // Get d-neighborhood of second node + canonically labeled graph
            cur2 = neweclass[it_classes[it2]][0]; // First element of neweclass
            get_neighborhood(sg, sub2, cur2, d);

            // Get canonical labelling
            t1 = clock();
            sparsenauty(&sub2, lab, ptn, orbits, &options, &stats, &cansub2);
            if(print_time_can_labelling)
               printf("Can time second: %d, %f\n", sub2.nv, ((double)(clock() - (double)t1))/CLOCKS_PER_SEC);
            
            // Compare on isomorphism and automorphism
            if(are_same_sg_can(&cansub1, &cansub2, cur1, cur2, v1_pos, lab, orbits)){
               neweclass[it_classes[it2]].push_back(it);
               class_id[it] = it_classes[it2]; // update class id of node it
               added = true;
               break;
            }
         }
      }
      // No class where node fits in, create a new class
      if(!added){
         class_id[it] = neweclass.size();
         neweclass.push_back({it});

         if(heuristic_choice == 1) class_key[std::make_pair(0, 0)].push_back(neweclass.size() - 1);
         else if(heuristic_choice == 2) degree_keys[degrees].push_back(neweclass.size() - 1); 
         else class_key[std::make_pair(n, m)].push_back(neweclass.size() - 1); // Default: n, m -> class nr
      }
   }

   // Free allocated memory
   SG_FREE(sub1);
   SG_FREE(sub2);
   SG_FREE(cansub1);
   SG_FREE(cansub2);
   
   if(print_statistics >= 3){
      printf("/final heuristic classes ");
      if(heuristic_choice == 2) printf("%ld\n", degree_keys.size());
      else printf("%ld\n", class_key.size());
      fflush(stdout);
   }
   return neweclass;
}

std::vector<int> twin_node_check(sparsegraph sg, std::map<int, std::set<int>> &twin_node_map){
   std::map<std::set<int>, int>targets; // Maps set of neighbours to a node
   std::set<int> neighbour_set;
   std::vector<int> node_set;

   // For each node, check if it has a twin node / clone
   for(size_t i = 0; i < sg.nv; i++){
         
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
         twin_node_count++;
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

std::vector< std::vector<int> >process_twin_nodes(std::vector< std::vector<int> >eq, std::map<int, std::set<int>>twin_node_map){
   std::vector< std::vector<int> >eq_new;
   std::vector<int> eq_class;
   std::set<int> nodes;

   for(auto it : eq){
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

std::vector< std::vector< int > > get_equivalence_classes(sparsegraph sg, const int d){
   std::vector< std::vector<int> >eq, eq_stat, temp, neweq;
   std::map<int, std::set<int>> twin_node_map;
   std::vector<int> node_set;
   size_t i, tel;
   size_t it_iso_checks = 0;
   size_t tot_iso_checks = 0;
   size_t it_can1 = 0;
   size_t tot_can1 = 0;
   clock_t t2;
   float t1;
   
   // Get set of all nodes: filter out twin nodes if any
   if(do_twin_node_check) 
      node_set = twin_node_check(sg, twin_node_map);
   else{
      for(i = 0; i < sg.nv; i++) 
         node_set.push_back(i);
   }
   eq.push_back(node_set);

   // Iteratively expand neighbourhood radius
   for(int i = 1; i <= d; i++){
      if(heuristic_choice == 0) i = d; // Naive computation -> start with i=d

      t2 = clock();
      tel = 0;
      it_iso_checks = 0;
      it_can1 = 0;

      // Iterate over equivalence classes and try to split them further
      for(auto it : eq){
         // Statistics
         iso_checks = 0;
         nr_can1 = 0;

         // Split current eq class
         temp = split_equivalence_class(sg, it, i);
         for(auto it2 : temp){
            neweq.push_back(it2);
         }

         // Print statistics for splitting eq class
         tel += 1;
         if(print_statistics >=3 && it.size() > 1){
            printf("/Finished class %ld / %ld\n", tel, eq.size()); 
            printf("/size %ld\n", it.size());
            printf("/time %f\n", ((double)(clock() - (double)t2))/CLOCKS_PER_SEC);
            printf("/can nbh1 %ld\n", nr_can1);
            printf("/iso checks %ld\n", iso_checks);
            printf("/new classes %ld\n", temp.size());
            fflush(stdout);
         }
         it_iso_checks += iso_checks;
         it_can1 += nr_can1;
      }
      // Update eq and clear neweq
      eq = neweq;
      neweq.clear();
      // Update statistics
      t1 += ((double)(clock() - (double)t2))/CLOCKS_PER_SEC;
      tot_iso_checks += it_iso_checks;
      tot_can1 += nr_can1;

      // Print statistics, get eqclass with twin nodes
      if(do_twin_node_check && (print_eq_class || print_statistics >=2))
         eq_stat = process_twin_nodes(eq, twin_node_map);
      else 
         eq_stat = eq;

      if(print_statistics >= 2){
         printf("N%d of %d:\n", i, d);
         printf("time it %d: %f\n", i, ((double)(clock() - (double)t2))/CLOCKS_PER_SEC);
         printf("k it %d: %d\n", i, get_k(eq_stat));
         printf("can nbh1 it %d: %ld\n", i, it_can1);
         printf("iso checks it %d: %ld\n", i, it_iso_checks);
         print_statistics_eq(eq_stat, i);
      }
      if(print_eq_class){
        print_equivalence_classes(eq_stat);
      }
      fflush(stdout);
      eq_stat.clear();
   }

   // Print statistics
   if(do_twin_node_check)
      eq = process_twin_nodes(eq, twin_node_map);
   if(print_statistics >= 1){
      printf("tot time: %f\n", t1);
      printf("final k: %d\n", get_k(eq));
      printf("tot can nbh1: %ld\n", tot_can1);
      printf("tot iso checks: %ld\n", tot_iso_checks);
      if(do_twin_node_check)
         printf("skipped twin nodes: %d\n", twin_node_count);      
      print_statistics_eq(eq, d);
   }
   fflush(stdout);
   return eq;
}

std::vector< std::vector< int > > split_equivalence_class_directed(sparsegraph sgo, sparsegraph sgi, std::vector <int> eclass, const int d){
   std::vector< std::vector< int > >neweclass = {};
   std::vector< int > it_classes;
   SG_DECL(sub1); SG_DECL(sub2);
   SG_DECL(cansub1); SG_DECL(cansub2);
   int cur1, cur2, n, m, nodes, edges, i, v1_pos;
   size_t tel;
   bool added;

   // Keys
   std::map< std::pair<int, int>, std::vector< int > > class_key; // Use nr. nodes + edges as map
   std::map< std::map<int, size_t>, std::vector< int > > outdeg_keys; // Use outdegree distribution as map
   std::map<int, size_t> outdegs;

   // Declare required arrays
   DYNALLSTAT(int, lab, lab_sz);
   DYNALLSTAT(int, ptn, ptn_sz);
   DYNALLSTAT(int, orbits, orbits_sz);
   DYNALLSTAT(int, map, map_sz);
   static DEFAULTOPTIONS_SPARSEGRAPH(options);
   options.getcanon = TRUE;
   statsblk stats;
   
   // Return if class has size 1, can not split further
   if(d == 0 || eclass.size() == 1){
      neweclass.push_back(eclass);
      return neweclass;
   }
   n = sgo.nv;
   m = sgo.elen;

   // Allocate space for neighborhood graphs
   SG_ALLOC(sub1, n, m, "malloc");
   SG_ALLOC(sub2, n, m, "malloc");

   tel = 0;
   
   // Iterate over nodes in current eclass, place them in new subclass
   for(auto it : eclass){
      added = false;
      cur1 = it;

      //Update:
      if(print_statistics >= 4 && (tel == 0 || tel % int((eclass.size() / 100) + 1) == 0)){
         printf("//node %ld / %ld\n", tel + 1, eclass.size());
         printf("//can nbh1 %ld\n", nr_can1);
         printf("//iso checks %ld\n", iso_checks);
         printf("//new classes %ld\n", neweclass.size());
         printf("//heuristic classes ");
         if(heuristic_choice == 2) printf("%ld\n", outdeg_keys.size());
         else printf("%ld\n", class_key.size());
         fflush(stdout);
      }
      tel += 1;

      // Get neighborhood of node to be placed
      get_neighborhood_directed(sgo, sgi, sub1, cur1, d);
      n = sub1.nv;
      m = sub1.nde;

      // Dynamic allocation for Traces
      DYNALLOC1(int, lab, lab_sz, n, "malloc");
      DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
      DYNALLOC1(int, orbits, orbits_sz,n, "malloc");
      DYNALLOC1(int, map, map_sz, n, "malloc");

      // Get keys using selected heuristic
      if(heuristic_choice == 0 || heuristic_choice == 1) // No heuristic
         it_classes = class_key[std::make_pair(0, 0)];
      else if(heuristic_choice == 2){ // Outdegree distribution
         get_degree_distribution_directed(sub1, outdegs);
         it_classes = outdeg_keys[outdegs];
      }
      else // Use nr. nodes and edges of neighborhood (default)
         it_classes = class_key[std::make_pair(n, m)];

      // Try to place node in one of the found classes
      if(it_classes.size() != 0){
         
         // Get canonically labeled graph and position of it in lab1 array
         sparsenauty(&sub1, lab, ptn, orbits, &options, &stats, &cansub1);
         for (i = 0; i < n; i++){
            if(lab[i] == cur1) {
               v1_pos = i; 
               break;
            }
         }
         nr_can1++;
         
         for(int it2 = 0; it2 < it_classes.size(); it2++){

            // Get d-neighborhood of second node + canonically labeled graph
            cur2 = neweclass[it_classes[it2]][0]; // A node in selected class of neweclass
            get_neighborhood_directed(sgo, sgi, sub2, cur2, d);
            sparsenauty(&sub2, lab, ptn, orbits, &options, &stats, &cansub2);
            
            // Compare on isomorphism and automorphism
            if(are_same_sg_can(&cansub1, &cansub2, cur1, cur2, v1_pos, lab, orbits)){
               neweclass[it_classes[it2]].push_back(it);
               added = true;
               break;
            }
         }
      }
      // No class where it fits in, create a new class
      if(!added){
         neweclass.push_back({it});
         if(heuristic_choice == 0 || heuristic_choice == 1) class_key[std::make_pair(0, 0)].push_back(neweclass.size() - 1);
         else if(heuristic_choice == 2) outdeg_keys[outdegs].push_back(neweclass.size() - 1); 
         else class_key[std::make_pair(n, m)].push_back(neweclass.size() - 1); // Default: n, m -> class nr
      }
   }

   if(print_statistics >= 3){
      printf("/final heuristic classes ");
      if(heuristic_choice == 2) printf("%ld\n", outdeg_keys.size());
      else printf("%ld\n", class_key.size());
      fflush(stdout);
   }

   // Free allocated memory
   SG_FREE(sub1);
   SG_FREE(sub2);
   SG_FREE(cansub1);
   SG_FREE(cansub2);
   
   return neweclass;
}

std::vector< std::vector< int > > get_equivalence_classes_directed(sparsegraph sgo, sparsegraph sgi, const int d){
   std::vector<int> all;
   std::vector< std::vector<int> >eq, temp, neweq;
   clock_t t2; 
   double t1 = 0.0;
   size_t i, n, tel;
   size_t it_iso_checks = 0;
   size_t tot_iso_checks = 0;
   size_t it_can1 = 0;
   size_t tot_can1 = 0;
   n = sgo.nv;

   // Start with class of all nodes (no information, all equivalent)
   for(int i = 0; i < n; i++){
      all.push_back(i);
   }
   eq.push_back(all);

   // Expand neighborhood
   for(int i = 1; i <= d; i++){
      t2 = clock();
      tel = 0;
      it_iso_checks = 0;
      it_can1 = 0;

      // Iterate over equivalence classes and try to split them further
      for(auto it : eq){
         iso_checks = 0;
         nr_can1 = 0;
         temp = split_equivalence_class_directed(sgo, sgi, it, i);
         for(auto it2 : temp){
            neweq.push_back(it2);
         }
         tel += 1;
         if(print_statistics >= 3 && it.size() > 1){
            printf("/Finished class %ld / %ld\n", tel, eq.size()); 
            printf("/size %ld\n", it.size());
            printf("/time %f\n", ((double)(clock() - (double)t2))/CLOCKS_PER_SEC);
            printf("/can nbh1 %ld\n", nr_can1);
            printf("/iso checks %ld\n", iso_checks);
            printf("/new classes %ld\n", temp.size());
            fflush(stdout);
         }
         it_iso_checks += iso_checks;
         it_can1 += nr_can1;
      }
      eq = neweq;
      neweq.clear();
      t1 += ((double)(clock() - (double)t2))/CLOCKS_PER_SEC;
      tot_iso_checks += it_iso_checks;
      tot_can1 += nr_can1;

      // Print statistics
      if(print_statistics >= 2){
         printf("N%d of %d:\n", i, d);
         printf("time it %d: %f\n", i, ((double)(clock() - (double)t2))/CLOCKS_PER_SEC);
         printf("k it %d: %d\n", i, get_k(eq));
         printf("can nbh1 it %d: %ld\n", i, it_can1);
         printf("iso checks it %d: %ld\n", i, it_iso_checks);
         print_statistics_eq(eq, i);
      }
      if(print_eq_class){
         print_equivalence_classes(eq);
      }
      fflush(stdout);
   }
   if(print_statistics >= 1){
      printf("tot time %f\n", t1);
      printf("tot can nbh1 %ld\n", tot_can1);
      printf("tot iso checks %ld\n", tot_iso_checks);
      if(print_statistics == 1){
         print_statistics_eq(eq, d);
      }
      fflush(stdout);
   }
   return eq;
}

int get_k(const std::vector< std::vector< int > > eclasses){
   int min_size;
   int cur_size;

   if(eclasses.size() == 0)
      return 0;
   
   min_size = eclasses[0].size();

   for(auto it : eclasses){
      cur_size = it.size();
      if(cur_size < min_size)
         min_size = cur_size;
   }
   return min_size;
}

void print_equivalence_classes(const std::vector< std::vector < int > > eclasses){
   printf("Start equivalence classes.\n");
   for(auto it : eclasses){
      for(auto it2 : it)
         printf("%d ", it2);
      printf("\n");
   }
   printf("End equivalence classes.\n");
}

void print_equivalence_classes_to_file(const std::vector < std::vector< int > > eclasses, char * file_name){
   FILE *fp;
   fp = fopen(file_name, "w+");
   
   for(auto it : eclasses){
      for(auto it2 : it)
         fprintf(fp, "%d ", it2);
      fprintf(fp, "\n");
   }
   fclose(fp);
}

void print_statistics_eq(const std::vector <std::vector <int > > eclasses, const int it){
   std::map<int, int> sizes;
   int max, min, cur;
   max = eclasses[0].size();
   min = eclasses[0].size();
   for(auto it : eclasses){
      cur = it.size();
      if(sizes.find(cur) == sizes.end()){
         sizes.emplace(cur, 0);
      }
      sizes[cur]++;
      if(cur < min) min = cur;
      if(cur > max) max = cur;
   }
   //printf("Equivalence partition statistics\n");
   //printf("It %d Min size: %d\n", it, min);
   //printf("It %d Max size: %d\n", it, max);
   //printf("All sizes (size, occurences):\n  ");
   printf("eq %d: ", it);
   for(auto it : sizes){
      printf("[%d, %d], ", it.first, it.second);
   }
   printf("\n");
   fflush(stdout);
}
