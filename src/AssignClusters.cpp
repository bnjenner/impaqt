#include <iostream>
#include <vector>

#include "utils.h"
#include "AnnotationList.h"
#include "ClusterList.h"
#include "AssignClusters.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Gene Assignment and Related Functions */

/////////////////////////////////////////////////////////////
/* Traversal Functions */

// Get Next Cloest Gene
GeneNode* get_closest_gene(const int &t, GeneNode *gene, ClusterNode *clust) {
	if (gene == NULL) { return NULL; }
	while (t > gene -> get_stop()) {
		gene = gene -> get_next();
		if (gene == NULL)  { break; }
		if (gene -> get_chrom() != clust -> get_contig_name()) {
			gene = NULL; break; // Not on same contig
		}
	}
	return gene;
}

/////////////////////////////////////////////////////////////
/* Assignment Functions */

// Resolve Read Assignment To Genes
void resolve_read_assignment(GeneNode *gene, const int &max, std::vector<size_t> &read_assignments) {

	if (gene == NULL && max != 0) {
		read_assignments[2] += 1; // Add to ambiguous

	} else if (max == 0) {
		read_assignments[1] += 1; // Add to Unassigned

	} else if (gene != NULL) {
		gene -> add_expression(1.0f); // add expression
		read_assignments[0] += 1;

 	} else {
 		std::cerr << "ERROR: No best overlap found for read assignment.\n";
 		throw "ERROR: No best overlap found for read assignment.";
 	}
}

// Resolve Read Assignment To Genes
void resolve_transcript_assignment(ClusterList *list, ClusterNode *node, GeneNode *gene, const int &max, const int &i) {
	
	long double expr = node -> get_transcript_expr(i);

	if (gene == NULL && max != 0) {
		node -> assign_ambiguous(i);
		list -> add_ambiguous_reads(expr);
		
	} else if (max == 0) {
		list -> add_unassigned_reads(expr);

	} else if (gene != NULL) {
		node -> assign_transcript(gene -> get_geneID(), i);
		gene -> add_expression(expr);
		list -> add_assigned_reads(expr);

	} else {
		std::cerr << "ERROR: No best overlap found for read assignment.\n"
				  << "Cluster: " << node -> get_contig_name() << ":" 
				  << node -> get_start() << "-" << node -> get_stop() << "\n";
 		throw "ERROR: No best overlap found for read assignment.";
 	}
}

/////////////////////////////////////////////////////////////
/* Overlapper Helper Functions */

int get_read_overlap(const int &a, const int &b, GeneNode *gene) {

	int match = 0;
	bool overlap_a, overlap_b;
	int m = gene -> get_exon_num();
	std::vector<int> exons = gene -> get_exon_vec();

	for (int i = 0; i < m; i++) {

		if (b < exons[(2*i)]) { break; }        // If exon past read
		if (exons[(2*i) + 1] < a) { continue; } // If read past exon

		overlap_a = check_point_overlap(a, exons[(2*i)], exons[(2*i)+1]);
		overlap_b = check_point_overlap(b, exons[(2*i)], exons[(2*i)+1]);

		// If perfect overlap
		if (overlap_a && overlap_b) {
			match = 2; break;
		
			// if partial overlap
		} else if (overlap_a || overlap_b) {
			match = 1; break;
		}
	}

	return match;
}

// Check the number of exons that overlap with transcript
int get_transcript_overlap(const std::vector<int> &transcript, GeneNode *gene) {

	int strand = gene -> get_strand();
	int matches = 0;
	int i = 0, j = 0;
	int n = transcript.size() / 2;
	int m = gene -> get_exon_num();
	std::vector<int> exons = gene -> get_exon_vec();
	bool overlap_a, overlap_b;

	// Iterate through clusters, see if match is complete 
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

					// prefer 3' most exon
					if (strand == 0 && j == m) {
						matches += 1;
					} else if (strand == 1 && j == 0) {
						matches += 1;
					}
				}
				j += 1;	
			}
		}
		i += 1;
	}

	return matches;
}

void compare_and_update_overlap(GeneNode *&gene, GeneNode *&best_gene, const int &overlap, int &max_overlap) {
	// If Better Overlap
	if (overlap > max_overlap) {
		max_overlap = overlap; best_gene = gene;

		 // If ambiguous
	} else if (overlap == max_overlap && max_overlap != 0) {
		best_gene = NULL;
	}
}

/////////////////////////////////////////////////////////////
/* Overlapper Functions */

// Assigning expression of transcripts to genes
void assign_transcripts_to_genes(ClusterNode *node, GeneNode *prev_gene, ClusterList *list,  AnnotationList &annotation, const int &t_num) {

	GeneNode *gene;
	GeneNode *best_gene;
	int t_start, t_stop, overlap, max_overlap;
	std::vector<std::vector<int>> transcripts = *(node -> get_transcripts());

	// Iterate through transcripts
	for (int i = 0; i < t_num; i++) {

		gene = prev_gene;
		best_gene = NULL;
		max_overlap = 0;
		
		// Get node start and stop
		t_start = transcripts[i][0];
		t_stop = transcripts[i][transcripts[i].size() - 1];

		// Check all possible genes
		while (t_stop >= gene -> get_start()) {
			
			overlap = get_transcript_overlap(transcripts[i], gene);
			compare_and_update_overlap(gene, best_gene, overlap, max_overlap);

			gene = gene -> get_next();
			if (gene == NULL || gene -> get_chrom() != node -> get_contig_name()) {
				break;
			}
		}

		resolve_transcript_assignment(list, node, best_gene, max_overlap, i);
	}
}

// Assigning expression of transcripts to genes
void assign_reads_to_genes(ClusterNode *node, GeneNode *prev_gene, ClusterList *list, AnnotationList &annotation) {

	GeneNode *gene;
	GeneNode *best_gene;
	int start, stop, overlap, max_overlap, index;

	int prev_read = -1;
	std::vector<size_t> read_assignments = {0, 0, 0}; // {Assigned, Unassigned, Ambiguous}, could probably make an array
	

	for (int i = 0; i < node -> get_vec_count(); i++) {

		index = node -> get_index_vec()[i];
		gene = prev_gene;

		// If onto new read, reset overlap stats
		if (index != prev_read) {

			if (prev_read != -1) {
				resolve_read_assignment(best_gene, max_overlap, read_assignments);
			}

			prev_read = index;
			max_overlap = 0;
			best_gene = NULL;
		}


		start = (node -> get_five_vec())[i];
		stop = (node -> get_three_vec())[i];

		// Check all possible genes
		while (stop >= gene -> get_start()) {
			overlap = get_read_overlap(start, stop, gene);
			compare_and_update_overlap(gene, best_gene, overlap, max_overlap);

			gene = gene -> get_next();
			if (gene == NULL || gene -> get_chrom() != node -> get_contig_name()) {
				break;
			}
		}
	}

	// Catch Last Read Assignment
	resolve_read_assignment(best_gene, max_overlap, read_assignments);

	// Increment Counts (necessary to do at once because of precision loss)
	for (int i = 0; i < 3; i++) {
		if (i == 0) {
			list -> add_assigned_singles(read_assignments[i]);
		} else if (i == 1) {
			list -> add_unassigned_singles(read_assignments[i]);
		} else {
			list -> add_ambiguous_singles(read_assignments[i]);
		}
	}
}

/////////////////////////////////////////////////////////////
/* Main Assignment Function */

void assign_to_genes(AnnotationList &annotation, ClusterList *list, const std::string &chrom, const int &strand) {
	
	int t_num, start;
	ClusterNode *node = list -> get_head(strand);
	GeneNode *prev_gene, *gene;

	if (ImpaqtArguments::Args.stranded == "reverse") {
		prev_gene = annotation.jump_to_chrom(chrom, !strand);
	} else {
		prev_gene = annotation.jump_to_chrom(chrom, strand);
	}

	// If no genes
	if (prev_gene == NULL) {
		list -> add_unassigned_singles(list -> get_passing_reads(strand));
		return;
	}

	while (node != NULL) {

		if (!(node -> is_skipped())) {

			t_num = node -> get_transcript_num();

			// If transcripts to assign
			if (t_num != 0) {
				start = node -> get_transcript_start();
				gene = get_closest_gene(start, prev_gene, node);
				if (gene != NULL) {
					assign_transcripts_to_genes(node, gene, list, annotation, t_num);
				}

			} else {

				node -> index_sort_vectors();
				start = (node -> get_five_vec())[0];
				gene = get_closest_gene(start, prev_gene, node);
				if (gene != NULL) {
					assign_reads_to_genes(node, gene, list, annotation);				
				}
				node -> empty_vectors();			}

			if (gene == NULL) {
				list -> add_unassigned_singles(node -> get_read_count());
			}
		}

		node = node -> get_next();
	}
}
