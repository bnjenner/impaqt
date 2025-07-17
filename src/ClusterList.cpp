#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <api/BamAux.h>
#include <api/BamReader.h>

#include "global_args.h"
#include "ClusterList.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Cluster Class Methods */

/////////////////////////////////////////////////////////////
/* Private Alignment Methods */

// Process CIGAR Strings
void ClusterList::calculate_splice(BamTools::BamAlignment &alignment, std::vector<int> &positions) {

	int n_offset = 0;
	int curr_pos = alignment.Position;
	int n = alignment.CigarData.size();
	positions.push_back(curr_pos);

	for (int i = 0; i < n; i++) {

		// If gapped alignment, get start and ends of neighbor aligned regions 
		if (alignment.CigarData[i].Type == 'N') {
			
			n_offset = alignment.CigarData[i].Length;
			curr_pos += n_offset;

		} else if (alignment.CigarData[i].Type == 'M') {
			curr_pos += alignment.CigarData[i].Length;

			// If gap detected and not last alignment
			if (n_offset != 0 && i != n - 1) {

				// If following D or I detected, skip, but only if no more gaps
				if ((alignment.CigarData[i + 1].Type == 'I' ||
					 alignment.CigarData[i + 1].Type == 'D') && 
					 i + 2 == n - 1) {
					continue;
				}

				positions.emplace_back(curr_pos - 1);
				positions.emplace_back(curr_pos - alignment.CigarData[i].Length);
				n_offset = 0;
			}
		}
		if (i == n - 1) { positions.emplace_back(curr_pos - 1); }
	}
}

// Check Read
bool ClusterList::read_check(const BamTools::BamAlignment &alignment) {

	if (alignment.IsDuplicate()) { return false; }
	if (!alignment.IsMapped()) { return false; }

	// Multimappers
	uint16_t NH_tag;
	alignment.GetTag("NH", NH_tag);
	if (!ImpaqtArguments::Args.nonunique_alignments && NH_tag > 1) {
		++ClusterList::multimapped_reads;
		return false;
	}

	// Exclude secondary alignment (do I need this?)
	if (!alignment.IsPrimaryAlignment() && !ImpaqtArguments::Args.nonunique_alignments) {
		++ClusterList::multimapped_reads;
		return false;
	}

	// If paired end, check propper pair
	if ((ImpaqtArguments::Args.library_type).compare("paired") == 0 && !alignment.IsProperPair()) {
		return false;
	}

	// Enfore MAPQ filter
	if (alignment.MapQuality < ImpaqtArguments::Args.mapq) {
		++ClusterList::low_quality_reads;
		return false;
	}

	return true;
}


/////////////////////////////////////////////////////////////
/* Private Node Methods */

// Create Empty Clusters
void ClusterList::initialize_strand(ClusterNode *&head, ClusterNode *&tail, const int strand, const int &zones) {
	int pos = 0;
	ClusterNode *node = new ClusterNode(pos,
	                                    strand,
	                                    ClusterList::window_size, 
	                                    ClusterList::contig_index, 
	                                    ClusterList::contig_name);
	head = node; tail = node;
	for (int i = 1; i < zones; i++) {
		pos += ClusterList::window_size;
		tail -> set_next(new ClusterNode(pos,
		                                 strand,
		                                 ClusterList::window_size, 
		                                 ClusterList::contig_index, 
		                                 ClusterList::contig_name));
		tail -> get_next() -> set_prev(tail);
		tail = tail -> get_next();
	}
}

// Delete Empty Nodes
void ClusterList::delete_nodes(ClusterNode *&c_node, ClusterNode *&t_head, ClusterNode *&t_tail) {

	ClusterNode *t_node;

	// If last node
	if (c_node -> get_next() == NULL) {
		t_node = c_node -> get_prev();

		// If curr_node is also not first node
		if (t_node != NULL) {
			t_node -> set_next(NULL);
		
		} else { t_head = NULL; }
		
		delete c_node;
		c_node = NULL; t_tail = t_node;
		return;
	}

	// If first node
	if (c_node -> get_prev() == NULL) {
		t_node = c_node -> get_next();

		t_node -> set_prev(NULL);
		delete c_node;

		c_node = t_node; t_head = c_node;
		return;
	}

	t_node = c_node -> get_prev();
	t_node -> set_next(c_node -> get_next());
	t_node -> get_next() -> set_prev(t_node);
	delete c_node;

	c_node = t_node -> get_next();
}

// Merge Neighboring Non-Zero Nodes
void ClusterList::merge_nodes(ClusterNode *&c_node, ClusterNode *&t_head, ClusterNode *&t_tail) {

	ClusterNode *n_node = c_node -> get_next();
	n_node -> shrink_vectors();

	ClusterNode *new_node = new ClusterNode(c_node, n_node);
	new_node -> set_prev(c_node -> get_prev());

	// If not first node
	if (c_node -> get_prev() != NULL) {
		c_node -> get_prev() -> set_next(new_node);
	
	} else { t_head = new_node; }

	// If not last node
	if (n_node -> get_next() != NULL) {
		new_node -> set_next(n_node -> get_next());
		new_node -> get_next() -> set_prev(new_node);
	
	} else { t_tail = new_node; }

	delete n_node; delete c_node;
	c_node = new_node;
}

// Delete every node in list
void ClusterList::delete_list() {
	ClusterNode *c_node = pos_head;
	ClusterNode *t_node = NULL;
	for (int i = 0; i < 2; i ++) {
		if (i != 0) { c_node = neg_head; t_node = NULL; }
		while (c_node != NULL) {
			t_node = c_node;
			c_node = c_node -> get_next();
			delete t_node;
		}
	}
}

/////////////////////////////////////////////////////////////
/* Cluster Functions */

// Create read clusters
bool ClusterList::create_clusters(BamTools::BamReader &inFile, BamTools::BamAlignment &alignment) {

	int t_5end, t_3end, t_strand;
	bool found_reads = false;
	std::vector<int> positions;
	ClusterNode *pos_node = ClusterList::get_head(0);
	ClusterNode *neg_node = ClusterList::get_head(1);
	ClusterNode *t_node = neg_node; // This needs to exist because BAM is ordered by left most position

	while (true) {

		if (!inFile.GetNextAlignment(alignment)) { break; }
		if (alignment.RefID > ClusterList::contig_index) { break; }

		total_reads += 1;

		if (ClusterList::read_check(alignment) == false) { continue; }

		found_reads = true;
		positions.clear();
		ClusterList::calculate_splice(alignment, positions); // Get Gapped Alignments

		// Process in Strand Specific way
		if (alignment.IsReverseStrand()) {

			t_5end = positions[positions.size() - 1];
			t_3end = positions[0];
			passing_neg_reads += 1;

			// Advance to correct node based on left position
			ClusterList::jump_to_cluster(neg_node, alignment.Position);
			t_node = neg_node;

			// Advance based on 5' position (right most)
			ClusterList::jump_to_cluster(t_node, t_5end);
			t_node -> add_alignment(positions);

		} else {

			t_5end = positions[0];
			t_3end = positions[positions.size() - 1];
			passing_pos_reads += 1;

			ClusterList::jump_to_cluster(pos_node, t_5end);
			pos_node -> add_alignment(positions);
		}
	}
	return found_reads;
}

// Combine clusters with nonzero neighbors
void ClusterList::collapse_clusters(int t_strand) {

	ClusterNode *t_head = ClusterList::get_head(t_strand);
	ClusterNode *t_tail = ClusterList::get_tail(t_strand);
	ClusterNode *node = t_head;

	while (node != NULL) {

		node -> shrink_vectors();

		// std::cerr << "Current: " << node -> get_start() << "-" << node -> get_stop() << ": " 
		// 		  << node -> get_read_count() << " reads\n";


		// Not Empty Node
		if (node -> get_read_count() != 0) {
			while (true) {

				// If no merging needed
				if (node -> get_next() == NULL) { break; }
				if (node -> get_next() -> get_read_count() == 0) { break; }

				// std::cerr << "\tMerging with: " << node -> get_next() -> get_start() 
				// 		  << " - " << node -> get_next() -> get_stop() << ": " 
				//  		  << node -> get_read_count() << " reads\n";"\n";


				ClusterList::merge_nodes(node, t_head, t_tail);
			}
			node = node -> get_next();

			// Empty
		} else { 
			// std::cerr << "\tDeleting\n";
			ClusterList::delete_nodes(node, t_head, t_tail);
		} // Also updates pointers

		// std::cerr << "\tHead: " << t_head -> get_start() << "-" << t_head -> get_stop() << ": "
		// 		  << t_head -> get_read_count() << " reads\n";
		// std::cerr << "\tTail: " << t_tail -> get_start() << "-" << t_tail -> get_stop() << ": "
		// 		  << t_tail -> get_read_count() << " reads\n";
	}

	// Adjust head and tail nodes by strand
	if (t_strand == 0) {
		ClusterList::pos_head = t_head; 
		ClusterList::pos_tail = t_tail;
	} else {
		ClusterList::neg_head = t_head; 
		ClusterList::neg_tail = t_tail;
	}
}


/////////////////////////////////////////////////////////////
/* Output Functions */

// Write Clusters to GTF File
void ClusterList::write_clusters_as_GTF(std::ofstream &gtfFile) {

	// Cancel if empty chromosome
	if (ClusterList::pos_head == NULL && ClusterList::neg_head == NULL) { return; }

	bool strand;
	ClusterNode *prev_pos = ClusterList::pos_head;
	ClusterNode *prev_neg = ClusterList::neg_head;
	ClusterNode *node = ClusterList::get_first_cluster(strand);

	// Iterate Through Clusters
	while (true) {

		// If no more clusters
		if (prev_pos == NULL && prev_neg == NULL) { break; }

		node -> ClusterNode::write_transcripts(gtfFile);

		// Get Next Cluster by position
		if (strand == 0) {
			node = ClusterList::get_next_cluster(node, prev_pos, prev_neg, strand);
		} else {
			node = ClusterList::get_next_cluster(node, prev_neg, prev_pos, strand);
		}
	}
}
