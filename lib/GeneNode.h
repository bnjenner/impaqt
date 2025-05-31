//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gene Node Class (node in a doubly linked list)
class GeneNode {

private:

	// Node Details
	std::string geneID;				   					// read ID of first read in cluster
	std::string chrom;					   				// chromosome index
	int strand = -1;					    			// standedness
	int start;											// beginning of window
	int stop;											// end of window
	int exons = 0;										// number of exons (or features)
	size_t read_count = 0;								// number of associated reads
	std::vector<int> exon_vec = {0, 0};					// vector for bounds

	// Links
	GeneNode *next = NULL;								// next ClusterNode
	GeneNode *prev = NULL;								// pevsious ClusterNode

	/////////////////////////////////////////////////////////////
	// If two exons overlap
	bool overlap(const int &e1, const int &e2, const int &t1, const int &t2) {
		if (e1 > t2 || e2 < t1) { return false; }
		return true;
	}

	/////////////////////////////////////////////////////////////
	// Append exon
	void append_exon(const int &t1, const int &t2) {
		int new_size = exon_vec.size() + 2;
		exon_vec.resize(new_size);
		exon_vec.at(new_size - 2) = t1;
		exon_vec.at(new_size - 1) = t2;
		exons += 1;
	}

	// Insert exon
	void insert_exon(const int &t1, const int &t2) {
		append_exon(t1, t2); // calls append (adds to end)
		std::sort(exon_vec.begin(), exon_vec.end()); // sort in ascending order
	}

	// Join exons
	void close_gap(const int &pos) {
        std::swap(exon_vec[(pos * 2) + 1], exon_vec.back()); exon_vec.pop_back(); // Remove end of current exon
        std::swap(exon_vec[(pos * 2) + 2], exon_vec.back()); exon_vec.pop_back(); // Remove beginning of next exon
        exons -= 1; // reduce number of exons
	}


public:

	/////////////////////////////////////////////////////////////
	// Empty
	GeneNode() {}

	// Initialize
	GeneNode(const std::string &t_geneID, const std::string &t_chrom, const std::string &t_strand, 
				const std::string &t_start, const std::string &t_stop) {
		// Set info
		geneID = t_geneID;
		chrom = t_chrom;
		strand = 0;
		if (t_strand == "-") { strand = 1; } // swtich if gene is on reverse

		// Set exons'
		exon_vec[0] = std::stoi(t_start) - 1;
		exon_vec[1] = std::stoi(t_stop) - 1;
		start = exon_vec[0]; 
		stop = exon_vec[1];
		exons += 1;
	}

	// Destroy
	~GeneNode() {
		next = NULL;
		prev = NULL;
	}

	/////////////////////////////////////////////////////////////
	// Get Details
	std::string get_chrom() { return chrom; }

	int get_strand() { return strand; }

	int get_start() { return start; }
	int get_stop() { return stop; }

	std::string get_geneID() { return geneID; }
	size_t get_read_count() { return read_count; }

	std::vector<int> get_exon_vec() { return exon_vec; }
	std::vector<int>* get_exon_ref() { return &exon_vec; }

	// Linking functions
	void set_next(GeneNode *node) { next = node; }
	void set_prev(GeneNode *node) { prev = node; }

	GeneNode* get_next() { return next; }
	GeneNode* get_prev() { return prev; }

	/////////////////////////////////////////////////////////////
	// Add exon to exon vector
	void add_region(const std::string &str_start, const std::string &str_stop) {

		const int n = exon_vec.size() / 2;
		const int t_start = std::stoi(str_start) - 1;
		const int t_stop = std::stoi(str_stop) - 1;

		
		// For all recorded exons
		for (int i = 0; i < n; i++) {

			// If intersection
			if (overlap(exon_vec[(2 * i)], exon_vec[(2 * i) + 1], t_start, t_stop)) {
				exon_vec[(2 * i)] = std::min(exon_vec[(2 * i)], t_start);
				exon_vec[(2 * i) + 1] = std::max(exon_vec[(2 * i) + 1], t_stop);

				// if gap now closed
				if (i != (n - 1) && exon_vec[(2 * i) + 1] >= exon_vec[(2 * i) + 2]) {
					close_gap(i);  // exon spans previous intron
				}
				break;
			}

			// If exon is further down stream
			if (t_start > exon_vec[(2 * i) + 1]) { 

				// if last exon
				if (i == (n - 1)) { 
					append_exon(t_start, t_stop);
				
				} else if (t_stop < exon_vec[(2 * i) + 2]) {
					insert_exon(t_start, t_stop);
				
				} else {
					continue;
				}
			}
		}
	}
};


