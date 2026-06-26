#include <iostream>
#include <vector>
#include <string>

#include "gtest/gtest.h"
#include "global_args.h"
#include "ClusterList.h"
#include "AnnotationList.h"
#include "AssignClusters.h"

// Globals (canonical GlobalArgs field order; see include/global_args.h)
ImpaqtArguments::GlobalArgs ImpaqtArguments::Args = {"",          // bam
                                                     "",          // index
                                                     "",          // annotation
                                                     1,           // threads
                                                     "forward",   // stranded
                                                     false,       // nonunique
                                                     1,           // mapq
                                                     2500,        // window size
                                                     1,           // min_count
                                                     25,          // count_percentage
                                                     50,          // epsilon
                                                     1.5,         // density threshold
                                                     false,       // isGFF
                                                     "exon",      // feature_tag
                                                     "UTR",       // utr_tag
                                                     "gene_id",   // feature_id
                                                     ""           // gtf_output
                                                    };


// Build a gene on chr1 (+) with three disjoint exons:
//   [99,199], [399,499], [699,799]  (0-based, inclusive; GeneNode subtracts 1)
static GeneNode make_three_exon_gene(const std::string &id) {
   GeneNode g(id, "chr1", "+", "100", "200");  // first exon -> [99,199]
   g.add_region("400", "500");                 // -> [399,499]
   g.add_region("700", "800");                 // -> [699,799]
   return g;
}


class AssignTest : public ::testing::Test {};


// get_read_overlap: 2 = both ends inside one exon, 1 = partial, 0 = miss.
TEST_F(AssignTest, ReadOverlapScoring) {
   GeneNode g = make_three_exon_gene("gene1");

   EXPECT_EQ(get_read_overlap(120, 180, &g), 2);  // fully inside exon 0
   EXPECT_EQ(get_read_overlap(420, 480, &g), 2);  // fully inside exon 1
   EXPECT_EQ(get_read_overlap(150, 250, &g), 1);  // 5' in exon 0, 3' past it
   EXPECT_EQ(get_read_overlap(150, 450, &g), 1);  // spans exon 0 -> exon 1 (still partial)
   EXPECT_EQ(get_read_overlap(250, 350, &g), 0);  // entirely in the gap
}


// get_transcript_overlap: sum of overlapping bases across exons.
TEST_F(AssignTest, TranscriptOverlapBases) {
   GeneNode g = make_three_exon_gene("gene1");

   // Segment wholly inside exon 0 -> 180 - 120 = 60 bases.
   std::vector<int> contained = {120, 180};
   EXPECT_EQ(get_transcript_overlap(contained, &g), 60);

   // Segment from inside exon 0, across the gap, into exon 1:
   //   exon 0 contributes 199 - 120 = 79; exon 1 fully covered = 499 - 399 = 100.
   std::vector<int> spanning = {120, 600};
   EXPECT_EQ(get_transcript_overlap(spanning, &g), 179);

   // Segment entirely in the gap -> no overlap.
   std::vector<int> miss = {220, 360};
   EXPECT_EQ(get_transcript_overlap(miss, &g), 0);
}


// compare_and_update_overlap: keep the strictly-better gene; a tie (non-zero)
// marks the result ambiguous (best_gene -> nullptr).
TEST_F(AssignTest, BestGeneSelectionAndAmbiguity) {
   GeneNode a("a", "chr1", "+", "100", "200");
   GeneNode b("b", "chr1", "+", "100", "200");
   GeneNode c("c", "chr1", "+", "100", "200");
   GeneNode *ga = &a, *gb = &b, *gc = &c;

   GeneNode *best = nullptr;
   int max_overlap = 0;

   compare_and_update_overlap(ga, best, 5, max_overlap);  // first real overlap
   EXPECT_EQ(best, ga);
   EXPECT_EQ(max_overlap, 5);

   compare_and_update_overlap(gb, best, 3, max_overlap);  // worse -> ignored
   EXPECT_EQ(best, ga);
   EXPECT_EQ(max_overlap, 5);

   compare_and_update_overlap(gc, best, 5, max_overlap);  // tie -> ambiguous
   EXPECT_EQ(best, nullptr);
   EXPECT_EQ(max_overlap, 5);
}


// get_closest_gene: walk forward until the gene's stop reaches t, stopping at
// the contig boundary.
TEST_F(AssignTest, ClosestGeneTraversal) {
   GeneNode a("a", "chr1", "+", "100", "200");   // stop = 199
   GeneNode b("b", "chr1", "+", "400", "500");   // stop = 499
   a.set_next(&b);  b.set_prev(&a);

   ClusterNode clust(100, 0, 2500, 0, "chr1");

   EXPECT_EQ(get_closest_gene(50,  &a, &clust), &a);       // before a -> a
   EXPECT_EQ(get_closest_gene(300, &a, &clust), &b);       // past a -> b
   EXPECT_EQ(get_closest_gene(600, &a, &clust), nullptr);  // past everything -> none

   // A gene on a different contig terminates the search.
   GeneNode other("z", "chr2", "+", "400", "500");
   a.set_next(&other);  other.set_prev(&a);
   EXPECT_EQ(get_closest_gene(300, &a, &clust), nullptr);
}


// resolve_read_assignment buckets a read into {assigned, unassigned, ambiguous}.
TEST_F(AssignTest, ResolveReadAssignmentBuckets) {
   GeneNode g("g", "chr1", "+", "100", "200");
   std::vector<size_t> ra = {0, 0, 0};

   resolve_read_assignment(&g, 2, ra);      // gene + overlap -> assigned
   EXPECT_EQ(ra[0], 1u);
   EXPECT_FLOAT_EQ(g.get_read_count(), 1.0f);

   resolve_read_assignment(nullptr, 2, ra); // overlap but no single best -> ambiguous
   EXPECT_EQ(ra[2], 1u);

   resolve_read_assignment(nullptr, 0, ra); // no overlap -> unassigned
   EXPECT_EQ(ra[1], 1u);

   resolve_read_assignment(&g, 0, ra);      // zero overlap dominates -> unassigned
   EXPECT_EQ(ra[1], 2u);
}


// End-to-end read assignment over a single gene: three reads -> two assigned,
// one unassigned (in the inter-exon gap), none ambiguous.
TEST_F(AssignTest, AssignReadsToGenesEndToEnd) {
   GeneNode g = make_three_exon_gene("gene1");
   g.set_next(nullptr);

   ClusterNode node(100, 0, 2500, 0, "chr1");
   node.add_alignment({120, 180}, {});  // read 0 -> inside exon 0  (assigned)
   node.add_alignment({420, 480}, {});  // read 1 -> inside exon 1  (assigned)
   node.add_alignment({250, 350}, {});  // read 2 -> in the gap     (unassigned)

   ClusterList *list = new ClusterList();  // heap-allocated; default dtor walks a
                                           // null list, so intentionally not freed
   assign_reads_to_genes(&node, &g, list);

   EXPECT_FLOAT_EQ((float)list->get_assigned_reads(),   2.0f);
   EXPECT_FLOAT_EQ((float)list->get_unassigned_reads(), 1.0f);
   EXPECT_FLOAT_EQ((float)list->get_ambiguous_reads(),  0.0f);
}
