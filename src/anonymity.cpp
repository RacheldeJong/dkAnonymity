// anonymity.cpp: Given a graph this code partitions the
// vertices based on d-k-anonymity.
// 
// Command line input: ./path-to/graph [optional arguments]
// Author: Rachel de Jong
// Date: 25-8-2021

#include "util.h"
#include "./graph/graphgen.h"
#include "./graph/graphutil.h"
#include "./measure/dk-anonymity.h"

// Read number of nodes from file
int read_n(const char* file_name){
   FILE* file;
   char c;
   int n = 0;

   if((file = fopen(file_name, "r")) == NULL){ 
      printf("Error: File %s could not be opened\n", file_name);
      return 0;
   }
   // Read n
   while (c != EOF){
      if(c == '\n' && (c=getc(file)) != '!') break;
      if(c == 'n' && (c=getc(file)) == '='){
         c = getc(file);
         while (ISDIGIT(c)){
            n *= 10;
            n += c - '0';
            c = getc(file);
         }
      }
      else c = getc(file);
   }
   return n;
}

// Print settings used
void print_info(sparsegraph sg1, const int directed, const int distance){
   printf("Info:\n");
   printf("- Graph contains %d nodes\n", sg1.nv);
   printf("- Graph contains %ld edges\n", sg1.nde);
   printf("- Directed: %d\n", directed);
   printf("- Distance: %d\n", distance);
   printf("- Heuristics choice = %d: ", heuristic_choice);
   
   if(heuristic_choice == 0) printf("Naive.\n");
   if(heuristic_choice == 1) printf("None.\n");
   else if(heuristic_choice == 2) printf("out-degree distribution\n");
   else printf("Nr. nodes and edges\n");
   
   printf("- Print equivalence classes = %d: \n", print_eq_class);
   printf("- Print statistics classes = %d: ", print_statistics);
   
   if(heuristic_choice <= 0) printf("none.\n\n");
   else if(print_statistics == 1) printf("Final statistics only\n\n");
   else if(print_statistics == 2) printf("Per iteration\n\n");
   else if(print_statistics == 3) printf("Per class split\n\n");
   else if(print_statistics >= 4) printf("All, including updates\n\n");
}

// Parse command line arguments
void parse_input(int argc, char* argv[], int & directed, int &distance){

   for(int i =1; i < argc; i++){
      if(strcmp(argv[i], "-dir") == 0){
         directed = atoi(argv[i + 1]);
         i++;
      }
      else if(strcmp(argv[i], "-d") == 0){
         distance = (atoi(argv[i + 1]));
         i++;
      }
      else if(strcmp(argv[i], "-h") == 0){
         heuristic_choice = (atoi(argv[i + 1]));
         i++;
      }
      else if(strcmp(argv[i], "-eq") == 0){
         print_eq_class = atoi(argv[i + 1]);
         i++;
      }
      else if(strcmp(argv[i], "-s") == 0){
         print_statistics = atoi(argv[i + 1]);
         i++;
      }
   }
}

// Main
int main(int argc, char *argv[]){
   SG_DECL(sg1); SG_DECL(sg2);
   int directed = 0;
   int distance = 5;
   int m, n = 0;
   char *input_file;
   heuristic_choice = 3;

   if(argc < 2){
      printf("Error not enough arguments: need input graph\n");
      exit(1);
   }
   
   input_file = argv[1];
   parse_input(argc, argv, directed, distance);
   printf("Graph file: %s\n", input_file);

   n = read_n(input_file);
   m = SETWORDSNEEDED(n);
   sg1 = read_graph_from_file(input_file, n);

   print_info(sg1, directed, distance);

   if(directed == 1){
     sg2 = get_ingoing_graph(sg1);
     get_equivalence_classes_directed(sg1, sg2, distance);
   }
   else{
     get_equivalence_classes(sg1, distance);
   }

   SG_FREE(sg1);
   SG_FREE(sg2);
   return 0;
}
