#include "gtest/gtest.h"
#include <iostream>
#include "api/BamReader.h"
#include "api/BamAux.h"
#include "parser.h"
#include "node.h"
#include "annotation.h"
#include "cluster.h"
#include "impaqt.h"
#include "queue.h"

class impactTest : public ::testing::Test {
public:

   // ImpactArguments args {
   //                       // Files
   //                       alignment_file;                     // sam or bam file
   //                       std::string index_file;             // index filed
   //                       std::string annotation_file;        // gff or gtf file

   //                       // Program
   //                       int threads;                        // threads

   //                       // Alignments
   //                       std::string library_type;           // library type (SE or PE)
   //                       std::string stranded;               // strandedness
   //                       bool nonunique_alignments;          // count primary and secondary alignments
   //                       int mapq;                           // minimum mapq score
   //                       int min_count;                      // min read count for dbscan
   //                       int count_percentage;                // read count percentage for core reads
   //                       int epsilon;                        // epsilon for dbscan

   //                       // Features
   //                       bool isGFF = false;                 // annotaiton file is gff
   //                       std::string feature_tag;            // name of feature tag
   //                       std::string feature_id;             // ID of feature

   //                       // Output
   //                       std::string gtf_output;             // name of output gtf file
   //                      };

};

TEST_F(impactTest, BasicCluster) { 
   std::cerr << "BasicCluster" << std::endl;
};
