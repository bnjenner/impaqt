//////////////////////////////////////
// General Node Class (node in a doubly linked list)
class Node {

public:

	// Node Variables
	int strand = -1;					    // standedness
	int clust_count = -1;				    // number of clusters
	int chrom_index;					    // chromosome number in index
	size_t read_count;					    // number of associated reads

	bool ishead;						    // is node head
	bool printed;						    // has node been printed

	// Clusters
	std::vector<int> count_vec{ -1 };       // vector of cluster read counts
	std::vector<int> clust_vec{ -1, -1 };   // vector of cluster start and stops (evens are starts, odds are ends)


	int get_start() { return clust_vec.at(0); }
	int get_stop() { return clust_vec.at(2 * (clust_count - 1) + 1); }

	// Check overlap
	bool check_overlap(int &t_start, int &t_stop, int &t_strand) {

		if (strand != t_strand) { return false; }

		for (int i = 0; i < clust_count; i++) {

			// check if beginning of read exists within a cluster
			if ((t_start >= clust_vec[i * 2]) && (t_start <= clust_vec[(i * 2) + 1])) {
				return true;

				// check if end of read exists within a cluster
			} else if ((t_stop >= clust_vec[i * 2]) && (t_stop <= clust_vec[(i * 2) + 1])) {
				return true;

				// in read spans cluster
			} else if ((t_start <= clust_vec[i * 2]) && (t_stop >= clust_vec[(i * 2) + 1])) {
				return true;
			}
		}

		return 0;
	}

	// insert gap
	void insert_read(int int_start, int int_stop, int new_pos, int new_val) {

		/*
			This function is a tiny bit cursed. Like it works, but
				also it's kind of the easiest way to do this. Insert new items
				and then sort them. Given the conditions met by the other function,
				it just works.
		*/

		// insert new start and stop
		clust_vec.reserve(clust_vec.size() + 2);
		clust_vec.emplace_back(int_start);
		clust_vec.emplace_back(int_stop);
		std::sort(clust_vec.begin(), clust_vec.end());

		// insert new count
		std::vector<int>::iterator it; // insert new position in count vector (probably not ideal)
		it = count_vec.begin() + new_pos;
		it = count_vec.insert(it, new_val);
	}

	// Delete gap
	void delete_gap(int t_index, int t_start, int t_stop) {

		int close_stop = 0;

		// Determine if clusters are joined by read
		for (int j = t_index + 1; j < clust_count; j++) {
			if (clust_vec[(j * 2)] > t_stop) { break; }
			close_stop = j;
		}

		if (close_stop != 0 && close_stop != t_index) {

			int start_index = (t_index * 2);
			int stop_index = (2 * close_stop) + 1;

			clust_vec.erase(clust_vec.begin() + start_index + 1, clust_vec.begin() + stop_index);

			int new_sum = 0;
			int offset = -1;
			std::vector<int> t_count_vec(clust_vec.size() / 2, -1);

			for (int x = 0; x < count_vec.size(); x++) {
				if (x < t_index || x > close_stop) {
					t_count_vec[x - std::max(offset, 0)] = count_vec[x];
				} else {
					new_sum += count_vec[x];
					offset++;
				}

			}

			t_count_vec[t_index] = new_sum;
			count_vec = t_count_vec;
		}
	}



	// check how read fits into clusters
	void modify_cluster(const int &t_start, const int &t_stop, const int t_count) {

		for (int i = 0; i < clust_count; i++) {

			// if read precedes cluster
			if (t_stop < clust_vec[(i * 2)]) {
				insert_read(t_start, t_stop, i, t_count);
				break;

				// if read follows cluster, go to next cluster... unless...
			} else if (t_start > clust_vec[(i * 2) + 1]) {

				if (i == clust_count - 1) { // if last cluster add
					insert_read(t_start, t_stop, i + 1, t_count);
					break;
				}

			} else {

				// Updates beginning of cluster to longest value betweeen end of cluster and end of read
				clust_vec[(i * 2)] = (clust_vec[(i * 2)] < t_start) ? clust_vec[(i * 2)] : t_start;

				// if end of read extends past first chunk
				if (t_stop > clust_vec[(i * 2) + 1]) {

					if (i != clust_count - 1) { delete_gap(i, t_start, t_stop); }

					// Updates end of cluster to longest value betweeen end of cluster and end of read
					clust_vec[(i * 2) + 1] = (clust_vec[(i * 2) + 1] > t_stop) ? clust_vec[(i * 2) + 1] : t_stop;
				}

				count_vec[i] += t_count;
				break;
			}
		}

		clust_count = clust_vec.size() / 2;
	}

};


class GeneNode : public Node {

public:

	// Links
	GeneNode *next = NULL;
	GeneNode *prev = NULL;

	std::string chrom = "";				   // name of chromosome
	std::string gene_id = "";			   // ID of gene (used in annotation graph)

	int overlap = 0;
	int max_overlap = 0;

	// Annotation Initialzed (for use in annotation)
	GeneNode(std::string temp_gene_id, int temp_strand, int begin, int end, std::string temp_chrom) {
		gene_id = temp_gene_id;
		strand = temp_strand;
		chrom = temp_chrom;
		read_count = 0;
		clust_count = 1;
		clust_vec = std::vector<int> {begin, end};
	}

	void set_next(GeneNode *node) { next = node; }
	void set_prev(GeneNode *node) { prev = node; }

	std::string get_chrom() { return chrom; }
	std::string get_gene() { return gene_id; }

	// Check genes and assigns overlap score (0 = no overlap, 1 = partial overlap, 2 = full overlap)
	int check_gene_overlap(int &temp_start, int &temp_stop, int &temp_strand) {

		if (strand != temp_strand) { return 0; } // is strand correct

		max_overlap = 0;

		for (int i = 0; i < clust_count; i++) {

			// if subcluster begins within exon
			if ((clust_vec[(i * 2)] <= temp_start) && (clust_vec[(i * 2) + 1] >= temp_start)) {

				// if subcluster is entirely within exon
				if (clust_vec[(i * 2) + 1] >= temp_stop) {
					max_overlap = 2;
				} else {
					max_overlap = std::max(max_overlap, 1);
				}

				// if subcluster starts before exon but overlaps with it
			} else if ((clust_vec[(i * 2)] > temp_start) && (clust_vec[(i * 2)] <= temp_stop)) {
				max_overlap = std::max(max_overlap, 1);
			}
		}
		return max_overlap;
	}
};


class ClusterNode : public Node {

public:

	// Links
	ClusterNode *next = NULL;
	ClusterNode *prev = NULL;

	std::vector<int> five_vec;			   // vector of read 5' site for transcript identification

	bool ambiguous = false;				   // is cluster ambigously assigned
	bool unassigned = true;			       // is cluster assigned
	std::string assigned_gene = "";		   // ID of gene assigned to cluster (used in cluster graph)

	// Empty
	ClusterNode() {}

	// Initialized (alignment)
	ClusterNode(BamTools::BamAlignment &alignment, int ref_num) {

		strand = alignment.IsReverseStrand();
		read_count = 1;
		chrom_index = ref_num;
		clust_vec[0] = alignment.Position;

		calculate_splice(alignment, clust_vec); // Calculate spliced alignments

		clust_count = clust_vec.size() / 2;
		count_vec[0] = 1;

		// push back cluster regions
		count_vec.reserve(count_vec.size() + clust_count - 1);
		for (int i = 1; i < clust_count; i++) {
			count_vec.emplace_back(1);
		}

		if (strand == 0) {
			five_vec.push_back(clust_vec[0]);
		} else {
			five_vec.push_back(clust_vec[1]);
		}

	}


	// Read cluster Initialized (region properties)
	ClusterNode(BamTools::BamAlignment &alignment, std::vector<int> &temp_vec, int ref_num) {
		strand = alignment.IsReverseStrand();
		read_count = 1;
		chrom_index = ref_num;
		clust_vec = temp_vec;
		clust_count = clust_vec.size() / 2;
		count_vec[0] = 1;

		// push back cluster regions
		count_vec.reserve(count_vec.size() + clust_count - 1);
		for (int i = 1; i < clust_count; i++) {
			count_vec.emplace_back(1);
		}


		if (strand == 0) {
			five_vec.push_back(clust_vec[0]);
		} else {
			five_vec.push_back(clust_vec[1]);
		}
	}

	// Transcript cluster Iniitialization
	ClusterNode(std::vector<int> &temp_vec, int ref_num, int temp_strand, int temp_count) {
		strand = temp_strand;
		read_count = temp_count;
		chrom_index = ref_num;
		clust_vec = temp_vec;
	}

	void set_next(ClusterNode *node) { next = node; }
	void set_prev(ClusterNode *node) { prev = node; }

	// Calculate splice
	void calculate_splice(BamTools::BamAlignment &alignment, std::vector<int> &temp_vec) {

		int pos = 1;
		int inc = 0;

		for (int i = 0; i < alignment.CigarData.size(); i++) {

			// If gap is encounterd, add splice, (may also need to check cigar string standards)
			if (alignment.CigarData[i].Type == 'N') {

				temp_vec[pos] = temp_vec[pos - 1] + inc - 1;
				temp_vec.reserve(temp_vec.size() + 2);

				// expand vector, slow, will improve (maybe)
				temp_vec.emplace_back(temp_vec[pos] + alignment.CigarData[i].Length + 1);
				temp_vec.emplace_back(-1);

				pos += 2;
				inc = 0;

				// If not gapped, add to start position
			} else if ((alignment.CigarData[i].Type != 'S') && (alignment.CigarData[i].Type != 'H') &&
			           (alignment.CigarData[i].Type != 'I')) {
				inc += alignment.CigarData[i].Length;
			}

			// if end of cigar string is reached
			if (i == alignment.CigarData.size() - 1) {
				temp_vec[pos] = temp_vec[pos - 1] + inc - 1;
			}
		}

	}


	// report cluster and counts
	int print(const std::string &contig_name, const ImpaqtArguments *parameters, int gene_count) {

		char s;
		std::string assignment;

		if ((clust_vec[0] == -1)) { return 0; }

		// Assign strand
		if ((parameters -> stranded).compare("forward") == 0) {
			s = (strand == 1) ? '-' : '+';
		} else {
			s = (strand == 1) ? '+' : '-';
		}

		// assignment
		if (unassigned) {
			assignment = "__unassigned";
		} else {
			if (ambiguous) {
				assignment = "__ambiguous";
			} else {
				assignment = "__assigned";
			}
		}

		std::ofstream outdata;
		outdata.open(parameters -> gtf_output, std::ios::out | std::ios::app);

		// Print "Gene" line, not contiguous
		outdata << contig_name << "\timpact\tcluster\t"
		        << clust_vec[0] + 1 << "\t" << clust_vec[((clust_count - 1) * 2) + 1] + 1
		        << "\t.\t" << s << "\t.\t"
		        << "gene_id \"impact." << contig_name << "." << gene_count << "\"; "
		        << "subclusters \"" << clust_count << "\"; "
		        << "counts \"" << read_count << "\"; "
		        << "assignment \"" << assignment << "\"; "
		        << "gene \"" << assigned_gene << "\";\n";

		// Iterate through clusters
		for (int i = 0; i < clust_count; i++) {
			// Print name, strand, and first start
			outdata << contig_name << "\timpact\tsubcluster\t" << clust_vec[(i * 2)] + 1 << "\t" << clust_vec[(i * 2) + 1] + 1;
			// Print the rest lol
			outdata << "\t.\t" << s << "\t.\t"
			        << "gene_id \"impact." << contig_name << "." << gene_count << "\"; "
			        << "subcluster_id \"impact." << contig_name << "." << gene_count << "." << i + 1 << "\"; "
			        << "counts \"" << count_vec[i] << "\";\n";
		}

		outdata.close();

		return 1;
	}
};