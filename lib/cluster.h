//////////////////////////////////////
// Cluster Class (really just a doubly linked list)
class ClusterList {

public:

	const ImpaqtArguments *parameters;	 	// parameters struct (found in parser.h)

	int chr_num = -1;						// chromosome number in index
	uint16_t NH_tag;						// NH tag to determine number of mappings
	std::string contig_name;				// name of contig

	// Links
	ClusterNode *head;
	ClusterNode *tail;
	ClusterNode temp;

	// Summary
	size_t multimapped_reads = 0;
	size_t total_reads = 0;

	// Empty
	ClusterList() {}

	// Initialize empty object
	void initialize(const int ref_num, const std::string ref_name, const ImpaqtArguments *args) {
		chr_num = ref_num;
		contig_name = ref_name;
		parameters = args;
	}

	// Check Read
	bool read_check(const BamTools::BamAlignment &alignment) {

			if (alignment.IsDuplicate()) { return false; }
			if (!alignment.IsMapped()) { return false; }

			alignment.GetTag("NH", NH_tag);
			if ((NH_tag > 1) && (!(parameters -> nonunique_alignments))) {
				multimapped_reads ++;
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

	// Set Head
	int set_head(BamTools::BamReader &inFile, BamTools::BamAlignment &alignment) {

		while (true) {

			if (!inFile.GetNextAlignment(alignment)) { return 0; }
			if (alignment.RefID > chr_num) { return 0; }

			total_reads ++;

			if (read_check(alignment) == false) { continue; }

			break;
		}

		// Add alignment to head node
		temp = ClusterNode(alignment, chr_num);
		temp.ishead = true;
		head = &temp;

		return 1;
	}

	// Create clusters of overlapping reads
	void create_clusters(BamTools::BamReader &inFile, BamTools::BamAlignment &alignment) {

		int temp_start;
		int temp_strand;
		size_t regions;
		std::vector<int> temp_vec = { -1, -1 };

		ClusterNode *new_node;
		ClusterNode *curr_node;
		tail = head;

		while (true) {

			curr_node = tail;

			if (!inFile.GetNextAlignment(alignment)) { break; }
			if (alignment.RefID > chr_num) { break; }

			total_reads ++;

			if (read_check(alignment) == false) { continue; } 

			temp_vec = {alignment.Position, -1};
			temp_start = temp_vec[0];
			temp_strand = alignment.IsReverseStrand();

			curr_node -> calculate_splice(alignment, temp_vec);
			regions = temp_vec.size() / 2;

			// check if alignment represents a new node (past first subcluster)
			if ((temp_vec[0] > curr_node -> clust_vec[1]) || (temp_strand != curr_node -> strand)) {

				new_node = new ClusterNode(alignment, temp_vec, chr_num);

				curr_node -> set_next(new_node);
				new_node -> set_prev(curr_node);
				curr_node = new_node;
				tail = curr_node;
				continue;
			}

			// find overlapping region
			while ((curr_node != NULL) && (temp_start != -1))  {

				for (int x = 0; x < regions; x++) {

					// Check if alignment overlaps with previous nodes
					if (curr_node -> check_overlap(temp_vec[(2 * x)], temp_vec[(2 * x) + 1], temp_strand)) {

						// add all clusters to vector
						for (int y = 0; y < regions; y++) {
							curr_node -> modify_cluster(temp_vec[(2 * y)], temp_vec[(2 * y) + 1], 1);
						}

						// if (regions < 2) {

						if (temp_strand == 0) {
							curr_node -> mid_vec.push_back(temp_vec[1]);
						} else {
							curr_node -> mid_vec.push_back(temp_vec[0]);
						}
							
							// curr_node -> mid_vec.push_back((temp_vec[0] + temp_vec[1]) / 2);
						// }
						
						curr_node -> read_count++;
						temp_start = -1;
						break;
					}
				}
				curr_node = curr_node -> prev;
			}
		}
		return;
	}

	// combines overlapping clusters in graph
	void collapse_clusters() {

		ClusterNode *curr_node = head;
		ClusterNode *next_node = NULL;
		ClusterNode *temp_node = NULL;
		ClusterNode *inter_node = NULL;

		// throw away variables, function only takes references :( will work on this
		int t_start;
		int t_stop;
		int t_strand;
		int t_overlap;
		int t_restart;

		while (curr_node != NULL) {

			next_node = curr_node -> next;
			temp_node = curr_node -> next;
			t_restart = 0;

			while (true) {

				if ((temp_node == NULL) ||
				        (temp_node -> get_start() > curr_node -> get_stop())) {
					if (t_restart == 0) {
						break;
					}

					t_restart = 0;
					temp_node = curr_node;

				} else {

					t_strand = temp_node -> strand;
					t_overlap = 0;

					for (int x = 0; x < temp_node -> clust_count; x++) {

						t_start = temp_node -> clust_vec[(x * 2)];
						t_stop = temp_node -> clust_vec[(x * 2) + 1];

						if (curr_node -> check_overlap(t_start, t_stop, t_strand) == 1) {
							t_overlap = 1;
							break;
						}

					}

					if (t_overlap == 1) {

						t_restart = 1;

						for (int x = 0; x < temp_node -> clust_count; x++) {
							t_start = temp_node -> clust_vec[(x * 2)];
							t_stop = temp_node -> clust_vec[(x * 2) + 1];
							curr_node -> modify_cluster(t_start, t_stop, temp_node -> count_vec[x]);
						}

						curr_node -> read_count += temp_node -> read_count;

						if (!(temp_node -> mid_vec.empty())) {
							curr_node -> mid_vec.reserve(curr_node -> mid_vec.size() + temp_node -> mid_vec.size());
							curr_node -> mid_vec.insert(std::end(curr_node -> mid_vec), 
														std::begin(temp_node -> mid_vec),
														std::end(temp_node -> mid_vec));
						}

						if ((temp_node -> next) != NULL) {
							(temp_node -> next) -> set_prev(temp_node -> prev);
						}

						(temp_node -> prev) -> set_next(temp_node -> next);
						inter_node = temp_node -> prev;

						if (temp_node == next_node) {
							next_node = temp_node -> next;
						}

						delete temp_node;
						temp_node = inter_node;
					}
				}

				temp_node = temp_node -> next;
			}

			curr_node = next_node;
		}
	}


	// print clusters
	void print_clusters() {

		int gene_count = 1;
		int printed;

		ClusterNode *curr_node = head;

		while (curr_node != NULL) {

			printed = curr_node -> print(contig_name, parameters, gene_count);
			curr_node -> printed = true;
			curr_node = curr_node -> next;

			if (printed == 1) { gene_count++; }
		}
	}
};