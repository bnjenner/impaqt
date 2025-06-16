#include "ClusterNode.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Cluster Class (really just a doubly linked list)
class ClusterList {

private:

	// List Details
	int chrom_index;					// chromosome number in index
	std::string contig_name;				// name of contig
	int chrom_length = 0;					// length of chromosome
	int window_size = 5000;					// window size

	// Links
	ClusterNode *pos_head;					// first positive ClusterNode
	ClusterNode *pos_tail;					// last positive ClusterNode
	ClusterNode *neg_head;					// first negative ClusterNode
	ClusterNode *neg_tail;					// last negative ClusterNode

	// Summary
	long double assigned_reads = 0.0;			// Assigned Transcript counts
	long double ambiguous_reads = 0.0;			// Unassigned Transcript counts
	long double unassigned_reads = 0.0;			// Ambigous Transcript counts
	size_t assigned_singles = 0;			// Assigned read counts
	size_t ambiguous_singles = 0;			// Unssigned read counts
	size_t unassigned_singles = 0;			// Ambigous read counts

	size_t multimapped_reads = 0;
	size_t low_quality_reads = 0;

	size_t total_reads = 0;
	size_t passing_pos_reads = 0;
	size_t passing_neg_reads = 0;
	size_t transcript_num = 0;



	/////////////////////////////////////////////////////////////
	// Calculate splice
	void calculate_splice(BamTools::BamAlignment &alignment, std::vector<int> &positions) {

		int n_offset = 0;
		int curr_pos = alignment.Position;
		int n = alignment.CigarData.size();
		positions.emplace_back(curr_pos);

		for (int i = 0; i < alignment.CigarData.size(); i++) {

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
	bool read_check(const BamTools::BamAlignment &alignment) {

		if (alignment.IsDuplicate()) { return false; }
		if (!alignment.IsMapped()) { return false; }

		// Multimappers
		uint16_t NH_tag;
		alignment.GetTag("NH", NH_tag);
		if ((!(ImpaqtArguments::Args.nonunique_alignments)) && (NH_tag > 1)) {
			++multimapped_reads;
			return false;
		}

		// Exclude secondary alignment (do I need this?)
		if (!alignment.IsPrimaryAlignment() && (!(ImpaqtArguments::Args.nonunique_alignments))) {
			++multimapped_reads;
			return false;
		}

		// If paired end, check propper pair
		if (((ImpaqtArguments::Args.library_type).compare("paired") == 0) && !alignment.IsProperPair()) {
			return false;
		}

		// Enfore MAPQ filter
		if (alignment.MapQuality < ImpaqtArguments::Args.mapq) {
			++low_quality_reads;
			return false;
		}

		return true;
	}


	void initialize_strand(ClusterNode *&head, ClusterNode *&tail, const int strand, const int &zones) {

		int pos = 0;
		ClusterNode *temp = new ClusterNode(pos,
							   				this -> window_size,
							   				strand, 
							   				this -> chrom_index,
							   				this -> contig_name);

		head = temp; tail = temp;
		for (int i = 1; i < zones; i++) {
			pos += this -> window_size;
			tail -> set_next(new ClusterNode(pos, 
											 this -> window_size, 
											 strand, 
											 this -> chrom_index,
											 this -> contig_name));
			tail -> get_next() -> set_prev(tail);
			tail = tail -> get_next();
		}
	}


	/////////////////////////////////////////////////////////////
	// Delete Empty Nodes (These node operations could likely simplified in functions...)
	void delete_nodes(ClusterNode *&curr_node, ClusterNode *&temp_head, ClusterNode *&temp_tail) {

		ClusterNode *temp_node;

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
	void merge_nodes(ClusterNode *&curr_node, ClusterNode *&temp_head, ClusterNode *&temp_tail) {

		ClusterNode *next_node = curr_node -> get_next();
		next_node -> shrink_vectors();

		ClusterNode *temp_node = new ClusterNode(*curr_node, *next_node);
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

	void delete_list() {
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
		pos_head = NULL;
		neg_head = NULL;
	}


public:

	/////////////////////////////////////////////////////////////
	// Empty
	ClusterList() {}

	// Destroy
	~ClusterList() { this -> delete_list(); }

	/////////////////////////////////////////////////////////////
	// Get Chrom Name
	std::string get_contig_name() { return contig_name; } // BNJ: 5/9/2025 - should probably keep name consistent with chrom, not contig

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
	long double get_assigned_reads() { return assigned_reads + (long double)assigned_singles; }
	long double get_unassigned_reads() { return unassigned_reads + (long double)unassigned_singles; }
	long double get_ambiguous_reads() { return ambiguous_reads + (long double)ambiguous_singles; }
	size_t get_multimapped_reads() { return multimapped_reads; }
	size_t get_low_quality_reads() { return low_quality_reads; }
	size_t get_total_reads() { return total_reads; }
	size_t get_pos_reads() { return passing_pos_reads; }
	size_t get_neg_reads() { return passing_neg_reads; }


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
		
	// Set Read Stats
	void add_assigned_reads(const long double &expr) { assigned_reads += expr; }
	void add_unassigned_reads(const long double &expr) { unassigned_reads += expr; }
	void add_ambiguous_reads(const long double &expr) { ambiguous_reads += expr; }
	void add_assigned_singles(const size_t &expr) { assigned_reads += expr; }
	void add_unassigned_singles(const size_t &expr) { unassigned_reads += expr; }
	void add_ambiguous_singles(const size_t &expr) { ambiguous_reads += expr; }

	/////////////////////////////////////////////////////////////
	// Initialize empty object
	void initialize(const int t_chrom_index, const std::string t_contig_name, const int t_chrom_length) {
		chrom_index = t_chrom_index;
		contig_name = t_contig_name;
		chrom_length = t_chrom_length;
		const int zones = (chrom_length / window_size) + 1;
		initialize_strand(pos_head, pos_tail, 0, zones); // Positive Strand
		initialize_strand(neg_head, neg_tail, 1, zones); // Negative Strand
	}

	// Find Nearest Region in List
	void jump_to_cluster(ClusterNode *&node, const int &pos) {
		while (node -> get_next() != NULL) {
			if (pos > node -> get_stop()) {
				node = node -> get_next();
			} else { break; }
		}	
	}

	/////////////////////////////////////////////////////////////
	// Create read clusters
	bool create_clusters(BamTools::BamReader &inFile, BamTools::BamAlignment &alignment) {

		int t_5end, t_3end, t_strand;
		bool found_reads = false;
		std::vector<int> positions;
		ClusterNode *pos_curr_node = get_head(0);
		ClusterNode *neg_curr_node = get_head(1);
		ClusterNode *neg_temp_node = neg_curr_node; // This needs to exist because the file is ordered according to the left most point

		while (true) {

			if (!inFile.GetNextAlignment(alignment)) { break; }
			if (alignment.RefID > chrom_index) { break; }

			total_reads += 1;

			if (read_check(alignment) == false) { continue; }

			found_reads = true;
			positions.clear();
			calculate_splice(alignment, positions); // Get Gapped Alignments

			// Process in Strand Specific way
			if (alignment.IsReverseStrand()) {

				t_5end = positions[positions.size() - 1];
				t_3end = positions[0];
				passing_neg_reads += 1;

				// Advance to correct node based on left position
				jump_to_cluster(neg_curr_node, alignment.Position);
				neg_temp_node = neg_curr_node;

				// Advance based on 5' position (right most)
				jump_to_cluster(neg_temp_node, t_5end);
				neg_temp_node -> add_alignment(positions);

			} else {

				t_5end = positions[0];
				t_3end = positions[positions.size() - 1];
				passing_pos_reads += 1;

				jump_to_cluster(pos_curr_node, t_5end);
				pos_curr_node -> add_alignment(positions);
			}
		}
		return found_reads;
	}

	/////////////////////////////////////////////////////////////
	// combines clusters with nonzero neighbors
	void collapse_clusters(int t_strand) {

		ClusterNode *temp_head = get_head(t_strand);
		ClusterNode *temp_tail = get_tail(t_strand);
		ClusterNode *curr_node = temp_head;

		while (curr_node != NULL) {

			curr_node -> shrink_vectors();

			// Not Empty Node
			if (curr_node -> get_read_count() != 0) {
				while (true) {

					// If no merging needed
					if (curr_node -> get_next() == NULL) { break; }
					if (curr_node -> get_next() -> get_read_count() == 0) { break; }

					merge_nodes(curr_node, temp_head, temp_tail);
				}
				curr_node = curr_node -> get_next();

				// Empty
			} else {
				delete_nodes(curr_node, temp_head, temp_tail);
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
	ClusterNode* get_first_cluster(ClusterNode *pos, ClusterNode *neg, bool &strand) {
		if (pos == NULL && neg != NULL) {
			strand = 1; return neg;
		} else if (neg == NULL && pos != NULL) {
			strand = 0; return pos;
		} else {
			if (pos -> get_start() < neg -> get_start()) {
				strand = 0; return pos;
			} else {
				strand = 1; return neg;
			}
		}
	}

	// Write Clusters to GTF File
	void write_clusters_as_GTF(std::ofstream &gtfFile) {

		// Cancel if empty chromosome
		if (pos_head == NULL && neg_head == NULL) { return; }

		bool strand;
		ClusterNode *prev_pos_node = pos_head;
		ClusterNode *prev_neg_node = neg_head;
		ClusterNode *curr_node = get_first_cluster(pos_head, neg_head, strand);

		// Iterate Through Clusters
		while (true) {

			// If no more clusters
			if (prev_pos_node == NULL && prev_neg_node == NULL) { break; }

			curr_node -> write_transcripts(gtfFile);

			// Strand switching conditions :(
			if (strand == 0) {
				prev_pos_node = curr_node -> get_next();

				// If positives exhausted, switch strands
				if (prev_pos_node == NULL) {
					curr_node = prev_neg_node; strand = 1;
					continue;
				}

				// If negatives exhausted, continue with positives
				if (prev_neg_node == NULL || prev_pos_node -> get_start() < prev_neg_node -> get_start()) {
					curr_node = prev_pos_node; strand = 0;

					// If positives are after negatives, switch strands
				} else { curr_node = prev_neg_node; strand = 1; }

			} else {
				prev_neg_node = curr_node -> get_next();

				// If negatives exhausted, switch strands
				if (prev_neg_node == NULL) { 
					curr_node = prev_pos_node; strand = 0;
					continue;
				}

				// If positives exhausted, continue with negatives
				if (prev_pos_node == NULL || prev_neg_node -> get_start() < prev_pos_node -> get_start()) {
					curr_node = prev_neg_node; strand = 1;

					// If negatives are after positives, switch strands
				} else { curr_node = prev_pos_node; strand = 0; }
			}
		}
	}

	/////////////////////////////////////////////////////////////
	// Functions for Test Suite
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
};