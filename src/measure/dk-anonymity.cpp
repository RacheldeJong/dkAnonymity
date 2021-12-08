#include "dk-anonymity.h"

int print_eq_class = 0;
int print_statistics = 1;
int heuristic_choice = 0;
size_t iso_checks, nr_can1;

bool are_same_sg(sparsegraph *sg1, sparsegraph *sg2, const int v1, const int v2){
   int i, n, m, target;
   n = sg1->nv;
   m = SETWORDSNEEDED(n);
   nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);

   // Declare required arrays
   DYNALLSTAT(int, lab1, lab1_sz);
   DYNALLSTAT(int, lab2, lab2_sz);
   DYNALLSTAT(int, ptn, ptn_sz);
   DYNALLSTAT(int, orbits, orbits_sz);
   DYNALLSTAT(int, map, map_sz);
   static DEFAULTOPTIONS_TRACES(options);
   TracesStats stats;

   // Check if number of nodes and edges is equal in graphs
   if(sg1->nv != sg2->nv || sg1->nde != sg2->nde){
      return false;
   }

   // Allocate space for arrays
   DYNALLOC1(int, lab1, lab1_sz, n, "malloc");
   DYNALLOC1(int, lab2, lab2_sz, n, "malloc");
   DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
   DYNALLOC1(int, orbits, orbits_sz,n, "malloc");
   DYNALLOC1(int, map, map_sz, n, "malloc");
   SG_DECL(cg1); SG_DECL(cg2);
   options.getcanon = TRUE;

   iso_checks++;

   // Nauty: Label sg1: results in cg1 and labelling in lab1; similarly sg2.
   // Also obtains orbits of each node in each graph
   Traces(sg1, lab1, ptn, orbits, &options, &stats, &cg1);
   Traces(sg2, lab2, ptn, orbits, &options, &stats, &cg2);

   // Compare canonically labelled graphs
   // If graphs are isomorph, then after Traces function sg1 == sg2 holds
   if (aresame_sg(&cg1,&cg2)){
      
      // Find node in sg2 isomorph with v1:
      // Each node i in the graph is mapped to a new label: lab1[i](G1) -> i (G'1)
      // lab1[i] =iso lab2[i] -> if v1 = lab1[i] then lab2[i] is a node in G2 isomorphic to v1 in G1
      for (i = 0; i < n; i++){
         if(lab1[i] == v1) {
            target = lab2[i]; 
            break;
         }
      }

      // Check if nodes are in the same orbit 
      // v1 =iso lab1[i] =iso lab2[i] and lab2[i] == target
      // -> v1 =iso target
      // if orbits[target] == orbits[v2] then target =iso v2 thus v1 =iso target =iso v2
      // -> v1 =iso v2
      if(orbits[v2] == orbits[target]){
         return true;
      }
   }
   DYNFREE(lab1, lab1_sz);
   DYNFREE(lab2, lab2_sz);
   DYNFREE(ptn, ptn_sz);
   DYNFREE(orbits, orbits_sz);
   DYNFREE(map, map_sz);

   return false;
}

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
   SG_DECL(sub1); SG_DECL(sub2);
   SG_DECL(cansub1); SG_DECL(cansub2);
   int cur1, cur2, n, m, nodes, edges, i, v1_pos;
   size_t tel;
   bool added;

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

      //Update:
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
            cur2 = neweclass[it_classes[it2]][0]; // First element of neweclass
            get_neighborhood(sg, sub2, cur2, d);
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

std::vector< std::vector< int > > get_equivalence_classes(sparsegraph sg, const int d){
   std::vector< std::vector<int> >eq, temp, neweq;
   std::vector<int> all;
   size_t i, n, tel;
   size_t it_iso_checks = 0;
   size_t tot_iso_checks = 0;
   size_t it_can1 = 0;
   size_t tot_can1 = 0;
   clock_t t2; 
   double t1 = 0.0;
   n = sg.nv;
   
   // Start with class of all nodes (no information, all equivalent)
   for(i = 0; i < n; i++){
      all.push_back(i);
   }
   eq.push_back(all);

   // Expand neighborhood
   for(int i = 1; i <= d; i++){
      if(heuristic_choice == 0) i = d; //
      t2 = clock();
      tel = 0;
      it_iso_checks = 0;
      it_can1 = 0;

      // Iterate over equivalence classes and try to split them further
      for(auto it : eq){
         iso_checks = 0;
         nr_can1 = 0;
         temp = split_equivalence_class(sg, it, i);
         for(auto it2 : temp){
            neweq.push_back(it2);
         }
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
      printf("tot time: %f\n", t1);
      printf("final k: %d\n", get_k(eq));
      printf("tot can nbh1: %ld\n", tot_can1);
      printf("tot iso checks: %ld\n", tot_iso_checks);
      if(print_statistics == 1)
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
