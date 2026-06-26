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
                                                     "../test/data/test.gtf",             // annotation
                                                     1,                                   // threads
                                                     "forward",                           // stranded
                                                     false,                               // nonunique
                                                     1,                                   // mapq
                                                     2500,                                // window size
                                                     25,                                  // min_count
                                                     5,                                   // count_percentage
                                                     150,                                 // epsilon
                                                     1.5,                                 // density threshold
                                                     false,                               // isGFF
                                                     "exon",                              // feature_tag
                                                     "UTR",                               // utr_tag
                                                     "gene_id",                           // feature_id
                                                     ""                                   // gtf_output
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

   static std::map<Path, int> paths;
   static std::vector<int> counts;
   static std::vector<std::vector<int>> transcripts;
   static std::vector<int> assign_vec_5, assign_vec_3;
   static std::map<int, std::vector<int>> regions_5, regions_3;

};


std::map<Path, int> impactTest::paths;
std::vector<int> impactTest::counts;
std::vector<std::vector<int>> impactTest::transcripts;
std::vector<int> impactTest::assign_vec_5; 
std::vector<int> impactTest::assign_vec_3;
std::map<int, std::vector<int>> impactTest::regions_5;
std::map<int, std::vector<int>> impactTest::regions_3;

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
   test_process -> create_clusters();
   test_process -> collapse_clusters();

   node = test_process -> get_clusters() -> get_head(0);
   expr = node -> get_read_count();
   points = node -> get_vec_count();
   min_counts = std::max((int)((float)expr * (((float)ImpaqtArguments::Args.count_percentage / 100.0))), 10);

   node -> point_sort_vectors();

   assign_vec_5 = dbscan(node, points, min_counts, regions_5, true);
   assign_vec_3 = dbscan(node, points, min_counts, regions_3, false);

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
   get_coordinates(paths,
                   regions_5, regions_3,
                   &transcripts, &counts);

   std::string result = "";
   for (const auto &p : transcripts) {
      for (const auto &pos : p) { result += std::to_string(pos) + ","; }
   }

   ASSERT_EQ(result, "4959707,4959824,4960962,4961018,4960962,4961094,4962137,4962291,4962137,4962284,");
};

// Test 2
TEST_F(impactTest, OverlapTranscripts) {

   overlap_clusters(node, transcripts, counts);
   report_transcripts(node, transcripts, counts);

   std::string result = "";
   for (const auto &p : transcripts) {
      for (const auto &pos : p) { result += std::to_string(pos) + ","; }
   }

   ASSERT_EQ(result, "4959707,4959824,4960962,4961094,4962137,4962291,");
};

// Test 3: regression for the >=10 cluster-index path bug.
//
// The old path key was a 2-char string: the producer did path = to_string(5') + '-'
// then path.at(1) = to_string(3')[0], and the consumer read substr(0,1)/substr(1,1).
// So a point linking 5' cluster 1 -> 3' cluster 10 was encoded "1-" then "11",
// colliding with a genuine 5'=1 -> 3'=1 link and silently dropping the real linkage.
// struct Path keys both as distinct integer pairs. This pins that fix.
TEST_F(impactTest, LargeClusterIndexPaths) {

   // A node whose only role is to report vec_count as the loop bound.
   ClusterNode big_node;
   std::vector<int> a5, a3;

   // 12 points linking 5' cluster 1 -> 3' cluster 1  (both single-digit)
   // 12 points linking 5' cluster 1 -> 3' cluster 10 (double-digit 3')
   // Under the old encoding both collapsed to "11"; they must stay distinct now.
   for (int i = 0; i < 12; i++) { a5.push_back(1); a3.push_back(1);  }
   for (int i = 0; i < 12; i++) { a5.push_back(1); a3.push_back(10); }
   big_node.update_vec_counts(a5.size());

   std::map<Path, int> local_paths;
   get_linked_clusters(&big_node, local_paths, a5, a3);

   // Two distinct surviving paths (each count 12 >= the threshold of 10),
   // not one merged "11" entry of 24.
   ASSERT_EQ(local_paths.size(), (size_t)2);
   ASSERT_EQ(local_paths[(Path{1, 1})],  12);
   ASSERT_EQ(local_paths[(Path{1, 10})], 12);

   std::map<int, std::vector<int>> r5, r3;
   r5[1]  = {1000, 1100};
   r3[1]  = {5000, 5100};
   r3[10] = {8000, 8100};   // would be read as r3[1] under the old substr(1,1) bug

   std::vector<std::vector<int>> local_transcripts;
   std::vector<int> local_counts;
   get_coordinates(local_paths, r5, r3, &local_transcripts, &local_counts);

   // Path{1,1} sorts before Path{1,10}; each keeps its own 3' region.
   ASSERT_EQ(local_transcripts.size(), (size_t)2);
   std::string result = "";
   for (const auto &p : local_transcripts) {
      for (const auto &pos : p) { result += std::to_string(pos) + ","; }
   }
   ASSERT_EQ(result, "1000,1100,5000,5100,1000,1100,8000,8100,");
};