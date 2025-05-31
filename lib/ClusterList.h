#include "ClusterNode.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Cluster Class (really just a doubly linked list)
class ClusterList {

private:

	// List Details
	int chrom_index;					    // chromosome number in index
	std::string contig_name;				// name of contig
	int chrom_length = 0;					// length of chromosome
	int window_size = 5000;					// window size

	// Links
	ClusterNode *pos_head;					// first positive ClusterNode
	ClusterNode *pos_tail;					// last positive ClusterNode
	ClusterNode *neg_head;					// first negative ClusterNode
	ClusterNode *neg_tail;					// last negative ClusterNode

	// Summary
	size_t total_reads = 0;
	size_t multimapped_reads = 0;
	size_t unassigned_reads = 0;
	size_t transcript_num = 0; 

	// For Checks
	uint16_t NH_tag;						// NH tag to determine number of mappings

	/////////////////////////////////////////////////////////////
	// Calculate splice
	int calculate_splice(BamTools::BamAlignment &alignment) {
		int end_pos = alignment.Position;
		for (int i = 0; i < alignment.CigarData.size(); i++) {
			// If not clipped and free of inserts, (may also need to check cigar string standards)
			if ((alignment.CigarData[i].Type != 'S') && (alignment.CigarData[i].Type != 'H') &&
			        (alignment.CigarData[i].Type != 'I')) {
				end_pos += alignment.CigarData[i].Length;
			}
		}
		return end_pos;
	}

	// Check Read
	bool read_check(const BamTools::BamAlignment &alignment) {

		if (alignment.IsDuplicate()) { return false; }
		if (!alignment.IsMapped()) { return false; }

		alignment.GetTag("NH", NH_tag);
		if ((NH_tag > 1) && (!(ImpaqtArguments::Args.nonunique_alignments))) {
			++multimapped_reads;
			return false;
		}

		// Exclude secondary alignment
		if (!alignment.IsPrimaryAlignment() && (!(ImpaqtArguments::Args.nonunique_alignments))) {
			return false;
		}

		// If paired end, check propper pair
		if (!alignment.IsProperPair() && ((ImpaqtArguments::Args.library_type).compare("paired") == 0)) {
			return false;
		}

		if (alignment.MapQuality <= ImpaqtArguments::Args.mapq) {
			return false;
		}

		return true;
	}

	/////////////////////////////////////////////////////////////
	// Delete Empty Nodes (These node operations could likely be cleaned up...)
	void delete_nodes(ClusterNode *&curr_node, ClusterNode *&temp_node,  ClusterNode *&temp_head, ClusterNode *&temp_tail) {
	
		// If last node
		if (curr_node -> get_next() == NULL) {
			temp_node = curr_node -> get_prev();

			// If curr_node is also not first node
			if (temp_node != NULL) {
				temp_node -> set_next(NULL);
			} else {
				temp_head = NULL;
			}

			temp_tail = temp_node;
			delete curr_node;
			curr_node = NULL;

			// If first node
		} else if (curr_node -> get_prev() == NULL) {
			temp_node = curr_node;
			curr_node = curr_node -> get_next();
			curr_node -> set_prev(NULL);
			temp_head = curr_node;
			delete temp_node;

			// else
		} else {
			temp_node = curr_node -> get_prev();
			temp_node -> set_next(curr_node -> get_next());
			temp_node -> get_next() -> set_prev(temp_node);

			delete curr_node;
			curr_node = temp_node -> get_next();
		}
	}

	// Merge Neighboring Non-Zero Nodes
	void merge_nodes(ClusterNode *&curr_node, ClusterNode *&temp_node, ClusterNode *&temp_head, ClusterNode *&temp_tail) {

		// Create new Node
		curr_node -> get_next() -> shrink_vectors();
		temp_node = new ClusterNode(curr_node -> get_start(),
		                            curr_node -> get_next() -> get_stop(),
		                            curr_node -> get_strand(),
		                            curr_node -> get_contig_name(),
		                            curr_node -> get_chrom_index(),
		                            curr_node -> get_read_count() + curr_node -> get_next() -> get_read_count(),
		                            curr_node -> get_five_vec(),
		                            curr_node -> get_next() -> get_five_vec(),
		                            curr_node -> get_three_vec(),
		                            curr_node -> get_next() -> get_three_vec());
		temp_node -> set_prev(curr_node -> get_prev());

		// If not first node
		if (curr_node -> get_prev() != NULL) {
			curr_node -> get_prev() -> set_next(temp_node);
		} else {
			temp_head = temp_node;
		}

		// If not last node
		if (curr_node -> get_next() -> get_next() != NULL) {
			temp_node -> set_next(curr_node -> get_next() -> get_next());
			temp_node -> get_next() -> set_prev(temp_node);
		} else {
			temp_tail = temp_node;
		}

		delete curr_node -> get_next();
		delete curr_node;
		curr_node = temp_node;
	}


public:

	/////////////////////////////////////////////////////////////
	// Empty
	ClusterList() {}

	// Destroy
	~ClusterList() {

		ClusterNode *curr_node = pos_head;
		ClusterNode *temp_node = NULL;

		while (curr_node != NULL) {
			temp_node = curr_node;
			curr_node = curr_node -> get_next();
			delete temp_node;
		}

		curr_node = neg_head;
		temp_node = NULL;

		while (curr_node != NULL) {
			temp_node = curr_node;
			curr_node = curr_node -> get_next();
			delete temp_node;
		}
	}

	/////////////////////////////////////////////////////////////
	// Get Chrom Name 
	//	BNJ: 5/9/2025 - should probably keep name consistent with chrom, not contig
	std::string get_contig_name() { return contig_name; }

	// Get Head Node
	ClusterNode* get_head(int t_strand) {
		if (t_strand == 0) { return pos_head; }
		return neg_head;
	}

	// Get Tail Node
	ClusterNode* get_tail(int t_strand) {
		if (t_strand == 0) { return pos_tail; }
		return neg_tail;
	}

	// Get Reads Stats
	size_t get_total_reads() { return total_reads; }
	size_t get_multimapped_reads() { return multimapped_reads; }
	
	// Get Transcript Number
	size_t get_transcript_num() {	
		ClusterNode *curr_node = pos_head;
		while (curr_node != NULL) {
			transcript_num += curr_node -> get_transcript_num();
			curr_node = curr_node -> get_next();
		}
		curr_node = neg_head;
		while (curr_node != NULL) {
			transcript_num += curr_node -> get_transcript_num();
			curr_node = curr_node -> get_next();
		}
		return transcript_num;
	}

	/////////////////////////////////////////////////////////////
	// Initialize empty object
	void initialize(const int t_chrom_index, const std::string t_contig_name, const int t_chrom_length) {

		// Set info
		chrom_index = t_chrom_index;
		contig_name = t_contig_name;
		chrom_length = t_chrom_length;

		ClusterNode *temp;
		int temp_pos = 0;
		int zones = (chrom_length / window_size) + 1; // Extend past length of chrom

		// Create Positive list
		temp = new ClusterNode(temp_pos, window_size, 0, chrom_index, contig_name);
		pos_head = temp; pos_tail = temp;
		for (int i = 1; i < zones; i++) {
			temp_pos += window_size;
			pos_tail -> set_next(new ClusterNode(temp_pos, window_size, 0, chrom_index, contig_name));
			pos_tail -> get_next() -> set_prev(pos_tail);
			pos_tail = pos_tail -> get_next();
		}

		// Create Negative list
		temp_pos = 0;
		temp = new ClusterNode(temp_pos, window_size, 1, chrom_index, contig_name);
		neg_head = temp; neg_tail = temp;
		for (int i = 1; i < zones; i++) {
			temp_pos += window_size;
			neg_tail -> set_next(new ClusterNode(temp_pos, window_size, 1, chrom_index, contig_name));
			neg_tail -> get_next() -> set_prev(neg_tail);
			neg_tail = neg_tail -> get_next();
		}
	}

	/////////////////////////////////////////////////////////////
	// Create read clusters
	bool create_clusters(BamTools::BamReader &inFile, BamTools::BamAlignment &alignment) {

		/*
		BNJ: 5/2/2025 - There's probably a better way to do this, will revisit...
		*/

		int t_5end, t_3end, t_strand;
		bool found_reads = false;
		ClusterNode *pos_curr_node = get_head(0);
		ClusterNode *neg_curr_node = get_head(1);
		ClusterNode *neg_temp_node = neg_curr_node; // This needs to exist because the file is ordered according to the left most point

		while (true) {

			if (!inFile.GetNextAlignment(alignment)) { break; }
			if (alignment.RefID > chrom_index) { break; }

			++total_reads;

			if (read_check(alignment) == false) { continue; }

			found_reads = true;

			if (alignment.IsReverseStrand()) {

				t_5end = calculate_splice(alignment);
				t_3end = alignment.Position;

				// Advance to correct node based on left position
				while (neg_curr_node -> get_next() != NULL) {
					if (alignment.Position >= neg_curr_node -> get_stop()) {
						neg_curr_node = neg_curr_node -> get_next();
					} else {
						break;
					}
				}

				neg_temp_node = neg_curr_node;

				// Get Correct node for 5' end (necessary because reverse strand)
				while (neg_temp_node -> get_next() != NULL) {
					if (t_5end >= neg_temp_node -> get_stop()) {
						neg_temp_node = neg_temp_node -> get_next();
					} else {
						break;
					}
				}

				// add read
				if (t_5end < neg_temp_node -> get_start()) {
					std::cerr << "ERROR: Alignment is not in a cluster. This should not happen.\n";
					throw "ERROR: Alignment is not in a cluster. This should not happen.";
				} else {
					neg_temp_node -> add_alignment(t_5end, t_3end);
				}

			} else {

				t_5end = alignment.Position;
				t_3end = calculate_splice(alignment);

				// Get Correct node
				while (pos_curr_node -> get_next() != NULL) {
					if (t_5end >= pos_curr_node -> get_stop()) {
						pos_curr_node = pos_curr_node -> get_next();
					} else {
						break;
					}
				}

				// Add read
				if (t_5end < pos_curr_node -> get_start()) {
					std::cerr << "ERROR: Alignment is not in a cluster. This should not happen.\n";
					throw "ERROR: Alignment is not in a cluster. This should not happen.";
				} else {
					pos_curr_node -> add_alignment(t_5end, t_3end);
				}
			}

		}

		return found_reads;
	}

	/////////////////////////////////////////////////////////////
	// combines clusters with nonzero neighbors
	void collapse_clusters(int t_strand) {

		ClusterNode *temp_node;
		ClusterNode *temp_head = NULL;
		ClusterNode *temp_tail = NULL;
		ClusterNode *curr_node = get_head(t_strand);

		while (curr_node != NULL) {

			curr_node -> shrink_vectors();

			// Not Empty Node
			if (curr_node -> get_read_count() != 0) {
				while (true) {

					// If no merging needed
					if (curr_node -> get_next() == NULL) { break; }
					if (curr_node -> get_next() -> get_read_count() == 0) { break; }

					merge_nodes(curr_node, temp_node, temp_head, temp_tail);
				}
				curr_node = curr_node -> get_next();

				// Empty
			} else {
				delete_nodes(curr_node, temp_node, temp_head, temp_tail);
			}
		}

		// Adjust head and tail nodes by strand
		if (t_strand == 0) {
			pos_head = temp_head; pos_tail = temp_tail;
		} else {
			neg_head = temp_head; neg_tail = temp_tail;
		}
	}

	/////////////////////////////////////////////////////////////
	// Print Clusters
	void print_clusters(int t_strand) {
		ClusterNode *curr_node = get_head(t_strand);
		while (curr_node != NULL) {
			std::cout << get_contig_name() << "\t"
			          << curr_node -> get_start() << "\t" << curr_node -> get_stop() << "\t"
			          << curr_node -> get_read_count() << "\n";
			curr_node = curr_node -> get_next();
		}
	}

	// Print clusters into strings for tests
	std::string string_clusters(int t_strand) {
		std::stringstream ss;
		ClusterNode *curr_node = get_head(t_strand);
		while (curr_node != NULL) {
			ss << get_contig_name() << "\t" 
			   << curr_node -> get_start() << "\t" << curr_node -> get_stop() << "\t" 
			   << curr_node -> get_read_count() << "\n";
			curr_node = curr_node -> get_next();
		}
		return ss.str();
	}

	// Write Clusters to GTF File
	void write_clusters_as_GTF(std::ofstream &gtfFile) {

		// Cancel if empty chromosome
		if (pos_head == NULL && neg_head == NULL) { return; }

		bool strand;
		ClusterNode *curr_node;
		ClusterNode *prev_pos_node = pos_head;
		ClusterNode *prev_neg_node = neg_head;

		// Get First Cluster
		if (pos_head == NULL && neg_head != NULL) {
			curr_node = neg_head; strand = 1;
		} else if (neg_head == NULL && pos_head != NULL) {
			curr_node = pos_head; strand = 0;
		} else {
			if (pos_head -> get_start() < neg_head -> get_start()) {
				curr_node = pos_head; strand = 0;
			} else {
				curr_node = neg_head; strand = 1;
			}
		}

		while (true) {

			// If no more clusters
			if (prev_pos_node == NULL && prev_neg_node == NULL) { break; }

			// Write transcripts into GTF
			curr_node -> write_transcripts(gtfFile);

			// Strand switching conditions :(
			//	BNJ:q 5/2/2025 - This is a bit messy, but it works, maybe put this into a function?
			if (strand == 0) {
				prev_pos_node = curr_node -> get_next();
				if (prev_pos_node != NULL) {
					if (prev_neg_node == NULL) {
						curr_node = prev_pos_node; strand = 0;
					} else if (prev_pos_node -> get_start() < prev_neg_node -> get_start()) {
						curr_node = prev_pos_node; strand = 0;
					} else {
						curr_node = prev_neg_node; strand = 1;
					}
				} else {
					curr_node = prev_neg_node; strand = 1;
				}

			} else {
				prev_neg_node = curr_node -> get_next();
				if (prev_neg_node != NULL) {
					if (prev_pos_node == NULL) {
						curr_node = prev_neg_node; strand = 1;
					} else if (prev_neg_node -> get_start() < prev_pos_node -> get_start()) {
						curr_node = prev_neg_node; strand = 1;
					} else {
						curr_node = prev_pos_node; strand = 0;
					}
				} else {
					curr_node = prev_pos_node; strand = 0;
				}
			}
		}
	}
};