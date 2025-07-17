#include <iostream>
#include <fstream>
#include <vector>
#include <thread>
#include <memory>
#include <condition_variable>
#include <chrono>

#include "gtest/gtest.h"
#include "global_args.h"
#include "impaqt.h"


// Globals
ImpaqtArguments::GlobalArgs ImpaqtArguments::Args = {"../test/data/test.bam",      // bam
                                                     "../test/data/test.bam.bai",  // index
                                                     "../test/data/test.gtf",      // annotation
                                                     1,                            // threads
                                                     "SE",                         // library type
                                                     "forward",                    // stranded
                                                     false,                        // nonunique
                                                     0,                            // mapq
                                                     5000,                         // window size
                                                     1,                            // min_count
                                                     25,                           // count_percentage
                                                     50,                           // epsilon
                                                     false,                        // isGFF
                                                     "exon",                       // feature_tag
                                                     "gene_id",                    // feature_id
                                                     ""                            // gtf_output
                                                    };


// Static Member Defintions
std::string Impaqt::alignment_file_name;
std::string Impaqt::index_file_name;
AnnotationList Impaqt::annotation;
std::unordered_map<int, std::string> Impaqt::contig_map;
std::unordered_map<int, int> Impaqt::contig_lengths;


// Test Class
class impactTest : public ::testing::Test {
public:

   Impaqt *test_process;

};

// Read In Results File as string
std::string read_test_file(std::string filename) {
   std::string file_contents;
   std::ifstream inputFile(filename);
   std::string line;
   while (std::getline(inputFile, line)) { file_contents += line + "\n"; }
   inputFile.close();
   return file_contents;
}

// Test 0
TEST_F(impactTest, ConstructClustList) {
   test_process = new Impaqt(0);
   test_process -> open_alignment_file();
   test_process -> set_chrom_order();
};


// Test 1
TEST_F(impactTest, SpliceTest) {

   ClusterList *cluster_list = new ClusterList();
   BamTools::BamReader SpliceFile;
   BamTools::BamAlignment alignment;

   // Open alignment file 
   if (!SpliceFile.Open("../test/data/SpliceTest.bam")) {
      std::cerr << "ERROR: Could not read alignment file: ../test/data/SpliceTest.bam\n";
      throw "ERROR: Make sure alignment file exists.";
   }

   int vec_count = 0;
   int read_count = 0;
   std::vector<int> positions;
   std::vector<int> five_vec;
   std::vector<int> three_vec;

   while (SpliceFile.GetNextAlignment(alignment)) {

      positions.clear();
      cluster_list -> calculate_splice(alignment, positions);

      // Add positions to 3 and 5 vec
      int n = positions.size();
      for (int i = 0; i < n; i++) {
         if (i % 2 == 0) {
            five_vec.push_back(positions[i]);
         } else {
            three_vec.push_back(positions[i]);
         }
      }

      ++read_count;
   }

   std::string answer = read_test_file("../test/data/SpliceTest.result");
   std::string result = "";
   for (int i = 0; i < five_vec.size(); i++) { result += std::to_string(five_vec[i]) + "," + std::to_string(three_vec[i]) + "\n"; }
  
   ASSERT_EQ(result, answer);
};

// Test 3
TEST_F(impactTest, BasicCluster) {
   test_process ->  create_clusters();
   std::string answer = read_test_file("../test/data/test_cluster.txt");
   ASSERT_EQ(test_process -> get_clusters() -> string_clusters(0), answer);
};

// Test 4
TEST_F(impactTest, BasicCollapse) {
   test_process -> collapse_clusters();
   std::string answer = read_test_file("../test/data/test_collapse.txt");
   ASSERT_EQ(test_process -> get_clusters() -> string_clusters(0), answer);
};

