#include "node.h"
#include "api/BamAux.h"
#include "api/BamReader.h"

//////////////////////////////////////
// Cluster Class (really just a doubly linked list)
class ClusterList {

private:

	// Paramters
	const ImpaqtArguments *parameters;	 	// parameters struct (found in parser.h)

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

	// For Checks
	bool initialized = false;

	// Summary
	size_t total_reads = 0;
	size_t multimapped_reads = 0;
	size_t unassigned_reads = 0;

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
		if ((NH_tag > 1) && (!(parameters -> nonunique_alignments))) {
			multimapped_reads++;
			return false;
		}

		// Exclude secondary alignment
		if (!alignment.IsPrimaryAlignment() && (!(parameters -> nonunique_alignments))) {
			return false;
		}

		// If paired end, check propper pair
		if (!alignment.IsProperPair() && ((parameters -> library_type).compare("paired") == 0)) {
			return false;
		}

		if (alignment.MapQuality < parameters -> mapq) {
			return false;
		}

		return true;
	}


public:

	// For Checks
	uint16_t NH_tag;						// NH tag to determine number of mappings

	// Temporary Variables
	int t_strand;
	int left_pos;
	int right_pos;

	// Empty
	ClusterList() {}

	// Destroy
	~ClusterList() {

		if (initialized) {

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
	}


	std::string get_contig_name() { return contig_name; }

	// Get Head Node
	ClusterNode* get_head(int t_strand) { 
		if (t_strand == 0) {
			return pos_head;
		} else {
			return neg_head;
		}
	}

	// Get Tail Node
	ClusterNode* get_tail(int t_strand) {
		if (t_strand == 0) {
			return pos_tail;
		} else {
			return neg_tail;
		}
	}

	// Get Reads Stats
	size_t get_total_reads() { return total_reads; }
	size_t get_multimapped_reads() { return multimapped_reads; }

	// Initialize empty object
	void initialize(const int t_chrom_index, const std::string t_contig_name, const int t_chrom_length, const ImpaqtArguments *args) {

		chrom_index = t_chrom_index;
		contig_name = t_contig_name;
		chrom_length = t_chrom_length;
		parameters = args;

		ClusterNode *temp;
		int temp_pos = 0;
		int zones = (chrom_length / window_size) + 1; // Extend past length of chrom

		// Create Positive list
		temp = new ClusterNode(temp_pos, window_size, 0, chrom_index);
		pos_head = temp;
		pos_tail = temp;
		for (int i = 1; i < zones; i++) {
			temp_pos += window_size;
			pos_tail -> set_next(new ClusterNode(temp_pos, window_size, 0, chrom_index));
			pos_tail -> get_next() -> set_prev(pos_tail);
			pos_tail = pos_tail -> get_next();
		}

		// Create Negative list
		temp_pos = 0;
		temp = new ClusterNode(temp_pos, window_size, 1, chrom_index);
		neg_head = temp;
		neg_tail = temp;
		for (int i = 1; i < zones; i++) {
			temp_pos += window_size;
			neg_tail -> set_next(new ClusterNode(temp_pos, window_size, 1, chrom_index));
			neg_tail -> get_next() -> set_prev(neg_tail);
			neg_tail = neg_tail -> get_next();
		}

		initialized = true;
	}

	// Create read clusters
	bool create_clusters(BamTools::BamReader &inFile, BamTools::BamAlignment &alignment) {

		int t_5end;
		int t_strand;
		bool found_reads = false;
		ClusterNode *pos_curr_node = get_head(0);
		ClusterNode *neg_curr_node = get_head(1);
		ClusterNode *neg_temp_node = neg_curr_node; // This needs to exist because the file is ordered according to the left most point

		while (true) {

			if (!inFile.GetNextAlignment(alignment)) { break; }
			if (alignment.RefID > chrom_index) { break; }

			total_reads ++;

			if (read_check(alignment) == false) { continue; }

			found_reads = true;

			if (alignment.IsReverseStrand()) {

				t_5end = calculate_splice(alignment);

				// Advance to correct node based on left position
				while (neg_curr_node -> get_next() != NULL) {
					if (alignment.Position >= neg_curr_node -> get_stop()) {
						neg_curr_node = neg_curr_node -> get_next();
					} else {
						break;
					}
				}

				neg_temp_node = neg_curr_node;

				// Get Correct node for 5' end
				while (neg_temp_node -> get_next() != NULL) {
					if (alignment.Position >= neg_temp_node -> get_stop()) {
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
					neg_temp_node -> add_alignment(t_5end);
				}

			} else {
				t_5end = alignment.Position;

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
					pos_curr_node -> add_alignment(t_5end);
				}
			}

		}

		return found_reads;
	}


	// combines clusters with nonzero neighbors
	void collapse_clusters(int t_strand) {

		ClusterNode *curr_node = get_head(t_strand);
		ClusterNode *temp_node;
		ClusterNode *temp_head = NULL;
		ClusterNode *temp_tail = NULL;

		while (curr_node != NULL) {
			
			// Empty Node
			if (curr_node -> get_read_count() == 0) {

				// If last node
				if (curr_node -> get_next() == NULL) {

					temp_node = curr_node -> get_prev();

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

			// Not empty
			} else {

				while (true) {

					if (curr_node -> get_next() == NULL) { break; }

					if (curr_node -> get_next() -> get_read_count() != 0) {

						temp_node = new ClusterNode(curr_node -> get_start(),
			                                        curr_node -> get_next() -> get_stop(),
			                                        curr_node -> get_strand(),
			                                        curr_node -> get_chrom_index(),
			                                        curr_node -> get_read_count() + curr_node -> get_next() -> get_read_count(),
			                                        curr_node -> get_five_vec(),
			                                        curr_node -> get_next() -> get_five_vec());

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

					} else {
						break;
					}

				}

				curr_node = curr_node -> get_next();
			}
		}

		if (t_strand == 0) {
			pos_head = temp_head;
			pos_tail = temp_tail;
		} else {
			neg_head = temp_head;
			neg_tail = temp_tail;
		}
	}


	// print clusters
	void print_clusters(int t_strand) {

		ClusterNode *curr_node = get_head(t_strand);

		while (curr_node != NULL) {

			std::cout << get_contig_name() << "\t"
			          << curr_node -> get_start() << "\t"
			          << curr_node -> get_stop() << "\t"
			          << curr_node -> get_read_count() << "\n";

			curr_node = curr_node -> get_next();
		}
	}
};