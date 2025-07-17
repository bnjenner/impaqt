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
ImpaqtArguments::GlobalArgs ImpaqtArguments::Args = {"../test/data/SpliceTest.bam",      // bam
                                                     "../test/data/SpliceTest.bam.bai",  // index
                                                     "../test/data/test.gtf",      // annotation
                                                     1,                            // threads
                                                     "SE",                         // library type
                                                     "forward",                    // stranded
                                                     false,                        // nonunique
                                                     1,                            // mapq
                                                     2500,                         // window size
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
class annotationTest : public ::testing::Test {
public:

   Impaqt *test_process;
   AnnotationList* annotation;

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
TEST_F(annotationTest, ConstructAnnoList) {
   test_process = new Impaqt(0);
   test_process -> add_annotation();
   annotation = test_process -> get_annotation();
};

// Test 1
TEST_F(annotationTest, ConstructGene) {
   std::string answer = test_process -> get_annotation() -> string_genes(0);
   ASSERT_EQ("chr1\tgene1\t24\t74\t99\t499\t574\t624\t649\t799\t\nchr3\tgene3\t34\t74\t99\t499\t574\t624\t649\t799\t\n", answer);
};

// Test 2
TEST_F(annotationTest, GetGenes) {
   bool strand;
   GeneNode *prev_pos = annotation -> get_head(0);
   GeneNode *prev_neg = annotation -> get_head(1);
   GeneNode *node = annotation -> get_first_gene(strand);
   std::string answer = node -> get_chrom() + ":" + node -> get_geneID() + ",";
   node = annotation -> get_next_gene(node, prev_pos, prev_neg, strand);
   answer += node -> get_chrom() + ":" + node -> get_geneID() + ",";
   node = annotation -> get_next_gene(node, prev_neg, prev_pos, strand);
   answer += node -> get_chrom() + ":" + node -> get_geneID();
   ASSERT_EQ("chr1:gene1,chr2:gene2,chr3:gene3", answer);
};


