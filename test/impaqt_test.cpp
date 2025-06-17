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
ImpaqtArguments::GlobalArgs ImpaqtArguments::Args = {"../test/data/test.bam",      // bam
                                                     "../test/data/test.bam.bai",  // index
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
std::string Impaqt::index;
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


// Test 1
TEST_F(impactTest, BasicAnnotation) {
   test_process -> add_annotation();
   std::string answer = test_process -> get_annotation() -> string_genes(0);
   ASSERT_EQ("chr1\tgene1\t49\t149\t174\t324\t374\t499\t524\t574\t599\t899\t\n", answer);
};

// // Test 2
// TEST_F(impactTest, BasicCluster) {
//    test_process ->  create_clusters();
//    std::string answer = read_test_file("../test/data/test_cluster.txt");
//    ASSERT_EQ(test_process -> get_clusters() -> string_clusters(0), answer);
// };

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
