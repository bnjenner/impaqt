#include "gtest/gtest.h"
#include <iostream>
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
                           0,                   // count_percentage
                           0,                   // epsilon
                           false,               // isGFF
                           "exon",              // feature_tag
                           "gene_id",           // feature_id
                           ""};                 // gtf_output

};

TEST_F(impactTest, BasicCluster) { 

   Impaqt test_process(&args, 0);
   test_process.open_alignment_file();
   test_process.set_chrom_order();
   test_process.find_clusters();
   test_process.cluster_list.print_clusters(0);
   
};
