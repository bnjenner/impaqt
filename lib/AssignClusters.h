//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gene Assignment and Related Functions

// Get Next Cloest Gene
GeneNode* get_closest_gene(const int &t, GeneNode *curr_gene) {
	while (t > curr_gene -> get_stop()) {
		curr_gene = curr_gene -> get_next();
		if (curr_gene == NULL) { break; }
	}
	return curr_gene;
}


// Check the number of exons that overlap with transcript
int get_overlap(const std::vector<int> &transcript, GeneNode *gene) {

	int i = 0;
	int j = 0;
	int n = transcript.size() / 2;
	int m = gene -> get_exon_num();
	std::vector<int> exons = gene -> get_exon_vec();

	int matches = 0;

	// Iterate through exons, see if match is complete 
	while (i < n) {

		j = 0;
		while (j < m) {

			// If Exon is past transcript 
			if (transcript[(2*i) + 1] < gene -> get_exon_vec()[(2*j)]) {
				break;

				// If transcript starts after exon ends, continue
			} else if (transcript[(2*i)] > gene -> get_exon_vec()[(2*j) + 1]) {
				j += 1;

				// If potential match
			} else {
				if (check_bounds(transcript[(2*i)], transcript[(2*i)+1], 
						gene -> get_exon_vec()[(2*j)], gene -> get_exon_vec()[(2*j)+1])) {
					matches += 1;
					break;
				}
				j += 1;	
			}
		}
		i += 1;
	}

	return matches;
}


// Main Assignment Function
void assign_transcripts_to_genes(AnnotationList &annotation, ClusterList &cluster_list, const std::string &chrom, const int &strand) {
	
	std::vector<std::vector<int>> transcripts;
	ClusterNode *curr_clust = cluster_list.get_head(strand);
	GeneNode *prev_gene = annotation.jump_to_chrom(chrom, strand);
	GeneNode *curr_gene = prev_gene;
	GeneNode *best_overlap;
	
	// If no genes
	if (curr_gene == NULL) {
		if (strand == 0) {
			float tmp_count = cluster_list.get_total_reads();
			tmp_count -= (float)cluster_list.get_multimapped_reads();
			tmp_count -= (float)cluster_list.get_low_quality_reads();
			cluster_list.add_unassigned_reads(tmp_count);
		}
		return;
	}

	int t_num;
	int t_start, t_stop;
	int overlap, max_overlap;
	float expr;

	while (curr_clust != NULL) {

		// Get Transcript Into
		transcripts = *(curr_clust -> get_transcripts());
		t_num = curr_clust -> get_transcript_num();

		// If transcripts to assign
		if (t_num != 0) {

			t_start = transcripts[0][0];

			// Advance to Potential Gene
			prev_gene = get_closest_gene(t_start, prev_gene);
			if (prev_gene == NULL) {
				cluster_list.add_unassigned_reads((float)curr_clust -> get_read_count());
				break;
			}

			// Iterate through transcripts
			for (int i = 0; i < t_num; i++) {

				max_overlap = -1;
				best_overlap = NULL;
				curr_gene = prev_gene;

				// Get cluster start and stop
				t_start = transcripts[i][0];
				t_stop = transcripts[i][transcripts[i].size() - 1];
				expr = curr_clust -> get_transcript_expr(i);

				// Check all possible genes
				while (t_stop >= curr_gene -> get_start()) {
					overlap = get_overlap(transcripts[i], curr_gene);

					// If Better Overlap
					if (overlap > max_overlap) {
						curr_clust -> assign_transcript(curr_gene -> get_geneID(), i);
						max_overlap = overlap;
						best_overlap = curr_gene;

						// If ambiguous
					} else if (overlap == max_overlap && max_overlap != 0) {
						best_overlap = NULL;
					}

					curr_gene = curr_gene -> get_next();
					if (curr_gene == NULL) { break; }
				}

				// If multiple assignments, mark ambiguous
				if (best_overlap == NULL && max_overlap != -1) {
					curr_clust -> assign_ambiguous(i);
					cluster_list.add_ambiguous_reads(expr);
					continue;
				}

				// If no overlap, add to unassigned reads
				if (max_overlap == -1) {
					cluster_list.add_unassigned_reads(expr);
					continue;
				}

				// Add Expression to best overlapping gene
				if (best_overlap != NULL) {
					best_overlap -> add_expression(expr);
					cluster_list.add_assigned_reads(expr);
				}
			}
		
		} else {
			cluster_list.add_unassigned_reads((float)(curr_clust -> get_read_count()));
		}

		curr_clust = curr_clust -> get_next();
	}

}