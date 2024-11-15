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
	ClusterNode temp; // This needs to exist for some reason... to be fixed

	// Checks
	bool hashead = false;

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

		hashead = true;

		return 1;
	}

	// Create clusters of overlapping reads
	void create_clusters(BamTools::BamReader &inFile, BamTools::BamAlignment &alignment) {

		int t_start;
		int t_strand;
		size_t t_regions;
		std::vector<int> t_vec = { -1, -1 };

		ClusterNode *new_node;
		ClusterNode *curr_node;
		tail = head;

		while (true) {

			curr_node = tail;

			if (!inFile.GetNextAlignment(alignment)) { break; }
			if (alignment.RefID > chr_num) { break; }

			total_reads ++;

			if (read_check(alignment) == false) { continue; }

			t_vec = {alignment.Position, -1};
			t_start = t_vec[0];
			t_strand = alignment.IsReverseStrand();

			curr_node -> calculate_splice(alignment, t_vec);
			t_regions = t_vec.size() / 2;

			// check if alignment represents a new node (past first subcluster)
			if ((t_vec[0] > curr_node -> clust_vec.at(1)) || (t_strand != curr_node -> strand)) {

				new_node = new ClusterNode(alignment, t_vec, chr_num);

				curr_node -> set_next(new_node);
				new_node -> set_prev(curr_node);
				curr_node = new_node;
				tail = curr_node;
				continue;
			}

			// find overlapping region
			while ((curr_node != NULL) && (t_start != -1))  {

				for (int x = 0; x < t_regions; x++) {

					// Check if alignment overlaps with previous nodes
					if (curr_node -> check_overlap(t_vec.at((2 * x)), t_vec.at((2 * x) + 1), t_strand)) {

						// add all clusters to vector
						for (int y = 0; y < t_regions; y++) {
							curr_node -> modify_cluster(t_vec.at((2 * y)), t_vec.at((2 * y) + 1), 1);
						}

						if (t_strand == 0) {
							curr_node -> five_vec.push_back(t_vec.at(0));
						} else {
							curr_node -> five_vec.push_back(t_vec.at(1));
						}


						curr_node -> read_count++;
						t_start = -1;
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

		/*
			The reason this function exists is that clusters can be created
				and then expaned upon later if they have massively gapped alignments.
				This means that they can end up overlapping with new clusters downstream.
				Potentially we stop tolerating this nonsense.
		*/

		ClusterNode *curr_node = head;
		ClusterNode *next_node = NULL;
		ClusterNode *temp_node = NULL;
		ClusterNode *inter_node = NULL;

		// throw away variables, function only takes references :( will work on this
		int t_start;
		int t_stop;
		int t_strand;

		bool t_overlap;
		bool t_restart;

		while (curr_node != NULL) {

			t_restart = false;

			next_node = curr_node -> next;
			temp_node = curr_node -> next;

			while (true) {

				// find possible overlapping node
				if ((temp_node == NULL) || (temp_node -> get_start() > curr_node -> get_stop())) {

					if (!t_restart) { break; }

					t_restart = false;
					temp_node = curr_node;

					// confirm overlap
				} else {

					t_strand = temp_node -> strand;
					t_overlap = false;

					for (int x = 0; x < temp_node -> clust_count; x++) {

						t_start = temp_node -> clust_vec[(x * 2)];
						t_stop = temp_node -> clust_vec[(x * 2) + 1];

						t_overlap = curr_node -> check_overlap(t_start, t_stop, t_strand);
						if (t_overlap) { break; }
					}

					if (t_overlap) {

						t_restart = true;

						// modify clusters
						for (int x = 0; x < temp_node -> clust_count; x++) {
							t_start = temp_node -> clust_vec[(x * 2)];
							t_stop = temp_node -> clust_vec[(x * 2) + 1];
							curr_node -> modify_cluster(t_start, t_stop, temp_node -> count_vec[x]);
						}

						// increase read count
						curr_node -> read_count += temp_node -> read_count;

						// add 5' ends to vector (order doesn't matter here cause of DBSCAN)
						if (!(temp_node -> five_vec.empty())) {
							curr_node -> five_vec.reserve(curr_node -> five_vec.size() + temp_node -> five_vec.size());
							curr_node -> five_vec.insert(std::end(curr_node -> five_vec),
							                             std::begin(temp_node -> five_vec),
							                             std::end(temp_node -> five_vec));
						}


						// properly link modified node
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

	///////////////////////////////////////////////////////////////////////////
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
	///////////////////////////////////////////////////////////////////////////

};