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
ImpaqtArguments::GlobalArgs ImpaqtArguments::Args = {"../test/data/dbscan_test.bam",      // bam
                                                     "../test/data/dbscan_test.bam.bai",  // index
                                                     "../test/data/test.gtf",      // annotation
                                                     1,                            // threads
                                                     "SE",                         // library type
                                                     "forward",                    // stranded
                                                     false,                        // nonunique
                                                     1,                            // mapq
                                                     2500,                         // window size
                                                     25,                           // min_count
                                                     5,                            // count_percentage
                                                     150,                          // epsilon
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
   int expr, points, min_counts;

   ClusterNode *node;

   static std::map<std::string, int> paths;
   static std::vector<int> counts;
   static std::vector<std::vector<int>> transcripts;
   static std::vector<int> assign_vec_5, assign_vec_3;
   static std::vector<std::vector<int>> assignments_5,  assignments_3;

};


std::map<std::string, int> impactTest::paths;
std::vector<int> impactTest::counts;
std::vector<std::vector<int>> impactTest::transcripts;
std::vector<int> impactTest::assign_vec_5; 
std::vector<int> impactTest::assign_vec_3;
std::vector<std::vector<int>> impactTest::assignments_5;
std::vector<std::vector<int>> impactTest::assignments_3;

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
TEST_F(impactTest, DBSCAN) {

   test_process = new Impaqt(0);
   test_process -> open_alignment_file();
   test_process -> set_chrom_order();
   test_process ->  create_clusters();
   test_process -> collapse_clusters();

   node = test_process -> get_clusters() -> get_head(1) -> get_next();
   expr = node -> get_read_count();
   points = node -> get_vec_count();
   min_counts = std::max((int)((float)expr * (((float)ImpaqtArguments::Args.count_percentage / 100.0))), 10);

   assign_vec_5 = dbscan(node, points, min_counts, assignments_5, true, false);
   assign_vec_3 = dbscan(node, points, min_counts, assignments_3, false, false);

   std::string result = "";
   for (const auto &p : assign_vec_5) { result += std::to_string(p); }
   for (const auto &p : assign_vec_3) { result += std::to_string(p); }
   result += "\n";

   // Generating the answer so I don't have to store this string;
   std::string answer = read_test_file("../test/data/dbscan.result");
   ASSERT_EQ(result, answer);
};

// Test 1
TEST_F(impactTest, GetCoordinates) {

   get_linked_clusters(node, paths, assign_vec_5, assign_vec_3);
   get_coordinates(node, paths,
                   assignments_5, assignments_3,
                   &transcripts, &counts);

   std::string result = "";
   for (const auto &p : transcripts) {
      for (const auto &pos : p) { result += std::to_string(pos) + ","; }
   }

   ASSERT_EQ(result, "61942684,61942962,61942684,61942878,61943136,61943440,61943093,61943440,");
};

// Test 2
TEST_F(impactTest, OverlapTranscripts) {

   overlap_clusters(node, transcripts, counts);
   report_transcripts(node, transcripts, counts);

   std::string result = ""; 
   for (const auto &p : transcripts) {
      for (const auto &pos : p) { result += std::to_string(pos) + ","; }
   }

   ASSERT_EQ(result, "61942684,61942962,61942684,61943440,");
};

// Test 3
TEST_F(impactTest, FullTest) {
   test_process = new Impaqt(0);
   test_process -> open_alignment_file();
   test_process -> set_chrom_order();
   test_process ->  create_clusters();
   test_process -> collapse_clusters();
   identify_transcripts(test_process -> get_clusters(), 1);

   std::string result = "";

   ClusterNode *node = test_process -> get_clusters() -> get_head(1);
   while (node != NULL) {
      for (const auto &transcript : *(node -> get_transcripts())) {
         for (const auto &pos : transcript) {
            std::cerr << pos << ",";
         }
         std::cerr << "\n";
      }
      node = node -> get_next();
   }

   ASSERT_EQ(result, "4871,5204,\n22376,22721,\n32383,32703,\n40154,40434,\n49903,50181,");
};