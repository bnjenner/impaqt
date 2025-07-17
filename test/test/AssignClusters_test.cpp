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
TEST_F(impactTest, Assignment) {};