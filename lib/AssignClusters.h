//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Gene Assignment and Related Functions */

/////////////////////////////////////////////////////////////
/* Traversal Functions */

// Get Next Cloest Gene
std::shared_ptr<GeneNode> get_closest_gene(const int &t, std::shared_ptr <GeneNode>gene, std::shared_ptr<ClusterNode> node);

/////////////////////////////////////////////////////////////
/* Assignment Functions */

// Resolve Read Assignment To Genes
void resolve_read_assignment(std::shared_ptr<GeneNode> gene, const int &max, std::vector<size_t> &read_assignments);

// Resolve Read Assignment To Genes
void resolve_transcript_assignment(std::shared_ptr<ClusterList> list, std::shared_ptr<ClusterNode> node, 
                                   std::shared_ptr<GeneNode> gene, const int &max, const int &i);

/////////////////////////////////////////////////////////////
/* Overlapper Helper Functions */

int get_read_overlap(const int &a, const int &b, std::shared_ptr<GeneNode> gene);

// Check the number of exons that overlap with transcript
int get_transcript_overlap(const std::vector<int> &transcript, std::shared_ptr<GeneNode> gene);

// Compare overlap of genes to best match so far
void compare_and_update_overlap(std::shared_ptr<GeneNode> &gene, std::shared_ptr<GeneNode> &best_gene, const int &overlap, int &max_overlap);

/////////////////////////////////////////////////////////////
/* Overlapper Functions */

// Assigning expression of transcripts to genes
void assign_transcripts_to_genes(std::shared_ptr<ClusterNode> node, std::shared_ptr<GeneNode> prev_gene, 
                                 std::shared_ptr<ClusterList> list, AnnotationList &annotation, const int &t_num);

// Assigning expression of transcripts to genes
void assign_reads_to_genes(std::shared_ptr<ClusterNode> node, std::shared_ptr<GeneNode> prev_gene, 
                           std::shared_ptr<ClusterList> list, AnnotationList &annotation);

/////////////////////////////////////////////////////////////
/* Main Assignment Function */
void assign_to_genes(AnnotationList &annotation, std::shared_ptr<ClusterList> list, const std::string &chrom, const int &strand);