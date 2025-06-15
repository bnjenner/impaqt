//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gene Assignment and Related Functions

// Get Next Cloest Gene
GeneNode* get_closest_gene(const int &t, GeneNode *curr_gene, ClusterNode *curr_clust) {
	if (curr_gene == NULL) { return NULL; }
	while (t > curr_gene -> get_stop()) {
		curr_gene = curr_gene -> get_next();
		if (curr_gene == NULL)  { break; }
		if (curr_gene -> get_chrom() != curr_clust -> get_contig_name()) {
			curr_gene = NULL; // Not on same contig
			break;
		}
	}
	return curr_gene;
}

// Resolve Read Assignment To Genes
void resolve_read_assignment(GeneNode *best_overlap, const int &max_overlap, std::vector<size_t> &read_assignments) {
	// If multiple assignments, mark ambiguous
	if (best_overlap == NULL && max_overlap != 0) {
		read_assignments[2] += 1; // Add to ambiguous

		// If no overlap, add to unassigned reads
    } else if (max_overlap == 0) {
		read_assignments[1] += 1; // Add to Unassigned

	// Add Expression to best overlapping gene
	} else if (best_overlap != NULL) {
		best_overlap -> add_expression(1.0f); // May be an issue
		read_assignments[0] += 1;
 	} else {
 		std::cerr << "ERROR: No best overlap found for read assignment.\n";
 		throw "ERROR: No best overlap found for read assignment.";
 	}
}


// Resolve Read Assignment To Genes
void resolve_transcript_assignment(ClusterList &cluster_list, ClusterNode *curr_clust, GeneNode *best_overlap, 
									const int &max_overlap, const int &i) {
	
	long double expr = curr_clust -> get_transcript_expr(i);

	// If multiple assignments, mark ambiguous
	if (best_overlap == NULL && max_overlap != 0) {
		curr_clust -> assign_ambiguous(i);
		cluster_list.add_ambiguous_reads(expr);

		// If no overlap, add to unassigned reads
	} else if (max_overlap == 0) {
		cluster_list.add_unassigned_reads(expr);

		// Add Expression to best overlapping gene
	} else if (best_overlap != NULL) {
		curr_clust -> assign_transcript(best_overlap -> get_geneID(), i);
		best_overlap -> add_expression(expr);
		cluster_list.add_assigned_reads(expr);
	} else {
 		std::cerr << "ERROR: No best overlap found for read assignment.\n"
 				  << "Cluster: " << curr_clust -> get_contig_name() << ":" 
 				  << curr_clust -> get_start() << "-" << curr_clust -> get_stop() << "\n";
 		throw "ERROR: No best overlap found for read assignment.";
 	}
}


int get_read_overlap(const int &a, const int &b, GeneNode *gene) {

	int match = 0;

	bool overlap_a, overlap_b;
	int m = gene -> get_exon_num();
	std::vector<int> exons = gene -> get_exon_vec();

	for (int i = 0; i < m; i++) {

		// If exon past read
		if (b < exons[(2*i)]) { break; }

		// If read past exon
		if (exons[(2*i) + 1] < a) { continue; }

		overlap_a = check_point_overlap(a, exons[(2*i)], exons[(2*i)+1]);
		overlap_b = check_point_overlap(b, exons[(2*i)], exons[(2*i)+1]);

		// If perfect overlap
		if (overlap_a && overlap_b) {
			match = 2;
			break;
		
			// if partial overlap
		} else if (overlap_a || overlap_b) {
			match = 1;
			break;
		}
	}

	return match;
}


// Check the number of exons that overlap with transcript
int get_transcript_overlap(const std::vector<int> &transcript, GeneNode *gene) {

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
			if (transcript[(2*i) + 1] < exons[(2*j)]) {
				break;

				// If transcript starts after exon ends, continue
			} else if (transcript[(2*i)] > exons[(2*j) + 1]) {
				j += 1;

				// If potential match
			} else {
				if (check_bounds(transcript[(2*i)], transcript[(2*i)+1], exons[(2*j)], exons[(2*j)+1])) {
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

// Assigning expression of transcripts to genes
void assign_transcripts_to_genes(ClusterNode *curr_clust, GeneNode *prev_gene, 
								 ClusterList &cluster_list, AnnotationList &annotation, const int &t_num) {

	GeneNode *curr_gene;
	GeneNode *best_overlap;
	int t_start, t_stop, overlap, max_overlap;
	std::vector<std::vector<int>> transcripts = *(curr_clust -> get_transcripts());


	// Iterate through transcripts
	for (int i = 0; i < t_num; i++) {

		max_overlap = 0;
		best_overlap = NULL;
		curr_gene = prev_gene;

		// Get cluster start and stop
		t_start = transcripts[i][0];
		t_stop = transcripts[i][transcripts[i].size() - 1];

		// Check all possible genes
		while (t_stop >= curr_gene -> get_start()) {

			overlap = get_transcript_overlap(transcripts[i], curr_gene);

			// If Better Overlap
			if (overlap > max_overlap) {
				max_overlap = overlap;
				best_overlap = curr_gene;

				// If ambiguous
			} else if (overlap == max_overlap && max_overlap != 0) {
				best_overlap = NULL;
			}

			curr_gene = curr_gene -> get_next();
			if (curr_gene == NULL || curr_gene -> get_chrom() != curr_clust -> get_contig_name()) {
				overlap = 0;
				break;
			}
		}

		resolve_transcript_assignment(cluster_list, curr_clust, best_overlap, max_overlap, i);
	}
}


// Assigning expression of transcripts to genes
void assign_reads_to_genes(ClusterNode *curr_clust, GeneNode *prev_gene,
						   ClusterList &cluster_list, AnnotationList &annotation) {

	GeneNode *curr_gene;
	GeneNode *best_overlap;
	int t_start, t_stop, overlap, max_overlap, index;

	std::vector<size_t> read_assignments = {0, 0, 0}; // {Assigned, Unassigned, Ambiguous}
	int prev_read = -1;

	for (int i = 0; i < curr_clust -> get_vec_count(); i++) {

		index = curr_clust -> get_index_vec()[i];
		curr_gene = prev_gene;

		// If onto new read, reset overlap stats
		if (index != prev_read) {

			if (prev_read != -1) {
				resolve_read_assignment(best_overlap, max_overlap, read_assignments);
			}

			prev_read = index;
			max_overlap = 0;
			best_overlap = NULL;
		}


		t_start = (curr_clust -> get_five_vec())[i];
		t_stop = (curr_clust -> get_three_vec())[i];

		// Check all possible genes
		while (t_stop >= curr_gene -> get_start()) {

			overlap = get_read_overlap(t_start, t_stop, curr_gene);

			// If Better Overlap
			if (overlap > max_overlap) {
				max_overlap = overlap;
				best_overlap = curr_gene;

				// If ambiguous
			} else if (overlap == max_overlap && max_overlap != 0) {
				best_overlap = NULL;
			}

			curr_gene = curr_gene -> get_next();
			if (curr_gene == NULL || curr_gene -> get_chrom() != curr_clust -> get_contig_name()) {
				overlap = 0;
				break;
			}
		}

	}

	// Catch Last Read Assignment
	resolve_read_assignment(best_overlap, max_overlap, read_assignments);

	// Increment Counts (necessary to do at once because of precision loss)
	for (int i = 0; i < 3; i++) {
		if (i == 0) {
			cluster_list.add_assigned_singles(read_assignments[i]);
		} else if (i == 1) {
			cluster_list.add_unassigned_singles(read_assignments[i]);
		} else {
			cluster_list.add_ambiguous_singles(read_assignments[i]);
		}
	}
}


// Main Assignment Function
void assign_to_genes(AnnotationList &annotation, ClusterList &cluster_list, const std::string &chrom, const int &strand) {
	
	ClusterNode *curr_clust = cluster_list.get_head(strand);
	GeneNode *prev_gene = annotation.jump_to_chrom(chrom, strand);

	// If no genes
	if (prev_gene == NULL) {

		if (strand == 0) { 
			cluster_list.add_unassigned_singles(cluster_list.get_pos_reads());
		} else {
			cluster_list.add_unassigned_singles(cluster_list.get_neg_reads());
		}
		return;
	}

	int t_num;
	int t_start;

	while (curr_clust != NULL) {

		t_num = curr_clust -> get_transcript_num();
		
		// If transcripts to assign
		if (t_num != 0) {

			t_start = (*(curr_clust -> get_transcripts()))[0][0];

			// Advance to Potential Gene
			prev_gene = get_closest_gene(t_start, prev_gene, curr_clust);
			if (prev_gene == NULL) {
				cluster_list.add_unassigned_singles(curr_clust -> get_read_count());
			} else {
				assign_transcripts_to_genes(curr_clust, prev_gene, cluster_list, annotation, t_num);
			}

		} else {

			curr_clust -> index_sort_vectors();

			t_start = (curr_clust -> get_five_vec())[0];

			// Advance to Potential Gene
			prev_gene = get_closest_gene(t_start, prev_gene, curr_clust);
			if (prev_gene == NULL) {
				cluster_list.add_unassigned_singles(curr_clust -> get_read_count());
			} else {
				assign_reads_to_genes(curr_clust, prev_gene, cluster_list, annotation);
			}
		}

		curr_clust = curr_clust -> get_next();
	}
}