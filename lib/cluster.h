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
	bool pos_hashead = false;
	bool neg_hashead = false;

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
	ClusterNode temp;

	// Empty
	ClusterList() {}

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

		int temp_pos = 0;
		int zones = (chrom_length / window_size) + 1; // Extent past length of chrom

		// Create Positive list
		temp = ClusterNode(temp_pos, window_size, 0, chrom_index);
		pos_head = &temp;
		pos_tail = &temp;
		for (int i = 1; i < zones; i++) {
			temp_pos += window_size;
			pos_tail -> set_next(new ClusterNode(temp_pos, window_size, 0, chrom_index));
			pos_tail -> get_next() -> set_prev(pos_tail);
			pos_tail = pos_tail -> get_next();
		}

		// Create Negative list
		temp_pos = 0;
		temp = ClusterNode(temp_pos, window_size, 1, chrom_index);
		neg_head = &temp;
		neg_tail = &temp;
		for (int i = 1; i < zones + 1; i++) {
			temp_pos += window_size;
			neg_tail -> set_next(new ClusterNode(temp_pos, window_size, 1, chrom_index));
			neg_tail -> get_next() -> set_prev(neg_tail);
			neg_tail = neg_tail -> get_next();
		}

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


	// combines overlapping clusters in graph
	void collapse_clusters() {}

	// 	/*
	// 		The reason this function exists is that clusters can be created
	// 			and then expaned upon later if they have massively gapped alignments.
	// 			This means that they can end up overlapping with new clusters downstream.
	// 			Potentially we stop tolerating this nonsense.
	// 	*/

	// 	ClusterNode *curr_node = head;
	// 	ClusterNode *next_node = NULL;
	// 	ClusterNode *temp_node = NULL;
	// 	ClusterNode *inter_node = NULL;

	// 	// throw away variables, function only takes references :( will work on this
	// 	int t_start;
	// 	int t_stop;
	// 	int t_strand;

	// 	bool t_overlap;
	// 	bool t_restart;

	// 	while (curr_node != NULL) {

	// 		t_restart = false;

	// 		next_node = curr_node -> next;
	// 		temp_node = curr_node -> next;

	// 		while (true) {

	// 			// find possible overlapping node
	// 			if ((temp_node == NULL) || (temp_node -> get_start() > curr_node -> get_stop())) {

	// 				if (!t_restart) { break; }

	// 				t_restart = false;
	// 				temp_node = curr_node;

	// 				// confirm overlap
	// 			} else {

	// 				t_strand = temp_node -> strand;
	// 				t_overlap = false;

					// for (int x = 0; x < temp_node -> clust_count; x++) {

					// 	t_start = temp_node -> clust_vec[(x * 2)];
					// 	t_stop = temp_node -> clust_vec[(x * 2) + 1];

					// 	t_overlap = curr_node -> check_overlap(t_start, t_stop, t_strand);
					// 	if (t_overlap) { break; }
					// }

					// if (t_overlap) {

					// 	t_restart = true;

					// 	// modify clusters
					// 	for (int x = 0; x < temp_node -> clust_count; x++) {
					// 		t_start = temp_node -> clust_vec[(x * 2)];
					// 		t_stop = temp_node -> clust_vec[(x * 2) + 1];
					// 		curr_node -> modify_cluster(t_start, t_stop, temp_node -> count_vec[x]);
					// 	}

					// 	// increase read count
					// 	curr_node -> read_count += temp_node -> read_count;

					// 	// add 5' ends to vector (order doesn't matter here cause of DBSCAN)
					// 	if (!(temp_node -> five_vec.empty())) {
					// 		curr_node -> five_vec.reserve(curr_node -> five_vec.size() + temp_node -> five_vec.size());
					// 		curr_node -> five_vec.insert(std::end(curr_node -> five_vec),
					// 		                             std::begin(temp_node -> five_vec),
					// 		                             std::end(temp_node -> five_vec));
					// 	}


					// 	// properly link modified node
					// 	if ((temp_node -> next) != NULL) {
					// 		(temp_node -> next) -> set_prev(temp_node -> prev);
					// 	}
					// 	(temp_node -> prev) -> set_next(temp_node -> next);

					// 	inter_node = temp_node -> prev;

					// 	if (temp_node == next_node) {
					// 		next_node = temp_node -> next;
	// 					}

	// 					delete temp_node;
	// 					temp_node = inter_node;
	// 				}
	// 			}

	// 			temp_node = temp_node -> next;
	// 		}

	// 		curr_node = next_node;
	// 	}
	// }

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