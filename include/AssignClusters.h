//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Gene Assignment and Related Functions */

/////////////////////////////////////////////////////////////
/* Traversal Functions */

// Get Next Cloest Gene
GeneNode* get_closest_gene(const int &t, GeneNode *gene, ClusterNode *clust);

/////////////////////////////////////////////////////////////
/* Assignment Functions */

// Resolve Read Assignment To Genes
void resolve_read_assignment(GeneNode *gene, const int &max, std::vector<size_t> &read_assignments);

// Resolve Read Assignment To Genes
void resolve_transcript_assignment(ClusterList *list, ClusterNode *cluster, GeneNode *gene, const int &max, const int &i);

/////////////////////////////////////////////////////////////
/* Overlapper Helper Functions */

int get_read_overlap(const int &a, const int &b, GeneNode *gene);

// Check the number of exons that overlap with transcript
int get_transcript_overlap(const std::vector<int> &transcript, GeneNode *gene);

// Compare overlap of genes to best match so far
void compare_and_update_overlap(GeneNode *&gene, GeneNode *&best_gene, const int &overlap, int &max_overlap);

/////////////////////////////////////////////////////////////
/* Overlapper Functions */

// Assigning expression of transcripts to genes
void assign_transcripts_to_genes(ClusterNode *cluster, GeneNode *prev_gene, ClusterList *list, 
								 AnnotationList &annotation, const int &t_num);

// Assigning expression of transcripts to genes
void assign_reads_to_genes(ClusterNode *cluster, GeneNode *prev_gene, ClusterList *list, AnnotationList &annotation);

/////////////////////////////////////////////////////////////
/* Main Assignment Function */
void assign_to_genes(AnnotationList &annotation, ClusterList *list, const std::string &chrom, const int &strand);