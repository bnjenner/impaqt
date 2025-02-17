//////////////////////////////////////
// Cluster Node Class (node in a doubly linked list)
class ClusterNode {

private:

	// Node Details
	int chrom_index;					    // chromosome number in index
	int strand = -1;					    // standedness
	int start;								// beginning of window
	int stop;								// end of window

	// Read Details
	std::string headID;				   		// read ID of first read in cluster
	size_t read_count;					    // number of associated reads
	std::vector<int> five_vec;				// vector for 5' ends

	// Links
	ClusterNode *next = NULL;				// next ClusterNode
	ClusterNode *prev = NULL;				// pevsious ClusterNode


public:

	// Empty
	ClusterNode() {}

	// Initialize
	ClusterNode(int &t_start, int &t_window_size, int t_strand, int &t_chrom_index) {

		start = t_start;
		stop = t_start + t_window_size;
		strand = t_strand;
		chrom_index = t_chrom_index;

		five_vec.resize(1000, 0);
		read_count = 0;

	}

	// Get Details
	void set_chrom_index(int t_chrom_index) { chrom_index = t_chrom_index; }
	int get_chrom_index() { return chrom_index; }

	bool isReverse() { return !strand; }

	int get_start() { return start; }
	int get_stop() { return stop; }

	std::string get_headID() { return headID; }
	size_t get_read_count() { return read_count; }
	std::vector<int> get_five_vec() { return five_vec; }

	// Linking functions
	void set_next(ClusterNode *node) { next = node; }
	void set_prev(ClusterNode *node) { prev = node; }
	ClusterNode* get_next() { return next; }
	ClusterNode* get_prev() { return prev; }

	void add_alignment(int t_5end) {
		five_vec.at(read_count) = t_5end;
		read_count += 1;

		if (read_count % 1000 == 0) {
			five_vec.resize(five_vec.size() + 1000, 0);
		}
	}
};

