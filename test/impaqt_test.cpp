#include "gtest/gtest.h"
#include <iostream>
#include <fstream>
#include "parser.h"
#include "annotation.h"
#include "impaqt.h"

class impactTest : public ::testing::Test {
public:

   ImpaqtArguments args = {"../test/data/test.bam",     // bam
                           "../test/data/test.bam.bai", // index
                           "../test/data/test.gtf",     // annotation
                           1,                   // threads
                           "SE",                // library type
                           "forward",           // stranded
                           false,               // nonunique
                           0,                   // mapq
                           1,                   // min_count
                           25,                  // count_percentage
                           50,                  // epsilon
                           false,               // isGFF
                           "exon",              // feature_tag
                           "gene_id",           // feature_id
                           ""};                 // gtf_output

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

TEST_F(impactTest, BasicCluster) { 

   test_process = new Impaqt(&args, 0);

   test_process -> open_alignment_file();
   test_process -> set_chrom_order();
   test_process -> find_clusters();

   std::string answer = read_test_file("../test/data/test_cluster.txt");
   ASSERT_EQ(test_process -> cluster_list.string_clusters(0), answer);
};

TEST_F(impactTest, BasicCollapse) { 

   test_process -> collapse_clusters();
   std::string answer = read_test_file("../test/data/test_collapse.txt");
   ASSERT_EQ(test_process -> cluster_list.string_clusters(0), answer);

};


TEST_F(impactTest, BasicDBSCAN) { 

   // test_process -> find_transcripts();
   // std::string answer = read_test_file("../test/data/test_collapse.txt");
   // ASSERT_EQ(test_process -> cluster_list.string_clusters(0), answer);

};
