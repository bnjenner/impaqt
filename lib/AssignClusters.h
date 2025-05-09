//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gene Assignment and Related Functions
void assign_transcripts_to_genes(AnnotationList &annotation, ClusterList &cluster_list, const int &strand) {

	/*
	Here is the procedure, iterate through clusters and genes
		If genes in cluster region
			If cluster has transcripts: assign to proper genes
			If cluster has no transcripts: assign individual 3' and 5' ends
				(if ends on same gene: unique, 
				if confliciting genes: ambiguous,
				if lone ends: unique,
				if no overlap: unassigned)

	*/

	// Get head and tail
	GeneNode *prev_gene = annotation.get_head(strand);
	GeneNode *curr_gene = prev_gene;
	ClusterNode *curr_clust = cluster_list.get_head(strand);

	while (curr_clust != NULL) {

		// // If not empty
		// if (curr_clust -> get_read_count() != 0) {
			
		// 	while (curr_gene != NULL) {

		// 		// If gene is in cluster
		// 		if (curr_clust -> get_start() >= curr_gene -> get_start() &&
		// 		    curr_clust -> get_stop() <= curr_gene -> get_stop()) {
		// 			curr_clust -> add_gene(curr_gene -> get_geneID());
		// 			break;
		// 		}
		// 		curr_gene = curr_gene -> get_next();
		// 	}
		// }

		// std::cout << curr_clust -> get_chrom_index() << ":"
		//           << curr_clust -> get_start() << "-" 
		//           << curr_clust -> get_stop() << "\t";
		// for (const auto &t : *(curr_clust -> get_transcripts())) {
		// 	for (const auto &region : t) {
		// 		std::cout << region << "\t";
		// 	}
		// 	std::cout << "\n";
		// }

		curr_clust = curr_clust -> get_next();
		// curr_gene = annotation.get_head(strand);
	}

}