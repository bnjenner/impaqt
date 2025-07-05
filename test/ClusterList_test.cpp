#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <string>
#include <vector>
#include <thread>

#include "gtest/gtest.h"
#include "api/BamAux.h"
#include "api/BamReader.h"
#include <global_args.h>
#include "ArgParser.h"
#include "impaqt.h"


// Globals
ImpaqtArguments::GlobalArgs ImpaqtArguments::Args = {"../test/data/SpliceTest.bam",      // bam
                                                     "../test/data/SpliceTest.bam.bai",  // index
                                                     "../test/data/test.gtf",      // annotation
                                                     1,                            // threads
                                                     "SE",                         // library type
                                                     "forward",                    // stranded
                                                     false,                        // nonunique
                                                     2500,                         // window size
                                                     0,                            // mapq
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

std::string read_test_file(std::string filename) {

   std::string file_contents;

   std::ifstream inputFile(filename);

   std::string line;
   while (std::getline(inputFile, line)) { file_contents += line + "\n"; }
   inputFile.close();

   return file_contents;
}

// Test 0
TEST_F(impactTest, OpenAlignmentFile) {
   test_process = new Impaqt(0);
   test_process -> open_alignment_file();
   test_process -> set_chrom_order();
};



// Test 2
TEST_F(impactTest, SpliceTest) {
   test_process ->  create_clusters();
   std::string answer = read_test_file("../test/data/SpliceTest.result");
   test_process -> get_clusters() -> get_head(0) -> shrink_vectors();
   std::vector<int> five_vec = test_process -> get_clusters() -> get_head(0) -> get_five_vec();
   std::vector<int> three_vec = test_process -> get_clusters() -> get_head(0) -> get_three_vec();

   std::string result = "";
   for (int i = 0; i < five_vec.size(); i++) { result += std::to_string(five_vec[i]) + "," + std::to_string(three_vec[i]) + "\n"; }

   ASSERT_EQ(result, answer);
};

// // Test 3
// TEST_F(impactTest, BasicCollapse) {
//    test_process -> collapse_clusters();
//    std::string answer = read_test_file("../test/data/test_collapse.txt");
//    ASSERT_EQ(test_process -> get_clusters() -> string_clusters(0), answer);
// };

// // Test 3
// TEST_F(impactTest, BasicDBSCAN) {

//    // test_process -> find_transcripts();
//    // std::string answer = read_test_file("../test/data/test_collapse.txt");
//    // ASSERT_EQ(test_process -> cluster_list.string_clusters(0), answer);

// };
