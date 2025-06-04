//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Cluster Node Class (node in a doubly linked list)
class ClusterNode {

private:

	// Node Details
	std::string contig_name;							// name of contig
	int chrom_index;									// chromosome number in index
	int strand = -1;									// standedness
	int start;											// beginning of window
	int stop;											// end of window

	// Read Details
	std::string headID;									// read ID of first read in cluster
	size_t read_count = 0;								// number of associated reads
	float total_core_points = 0;						// number of total core points
	std::vector<int> five_vec;							// vector for 5' ends
	std::vector<int> three_vec;							// vector for 3' ends

	// Links
	ClusterNode *next = NULL;							// next ClusterNode
	ClusterNode *prev = NULL;							// pevsious ClusterNode

	// Transcript Results
	int transcript_num = 0;								// number of transcripts identified
	std::vector<std::vector<int>> transcript_vec;		// vector of transcript regions
	std::vector<float> transcript_expression;			// vector of transcript expression
	std::vector<std::string> transcript_assignments;	// vector of transcript assignments


public:

	/////////////////////////////////////////////////////////////
	// Empty
	ClusterNode() {}

	// Initialize
	ClusterNode(int &t_start, int &t_window_size, int t_strand, int &t_chrom_index, std::string &t_contig_name) {
		start = t_start;
		stop = t_start + t_window_size;
		strand = t_strand;
		chrom_index = t_chrom_index;
		contig_name = t_contig_name;
		five_vec.resize(1000, 0);
		three_vec.resize(1000, 0);
	}

	// For combined Nodes
	ClusterNode(int t_start, int t_stop, int t_strand,
	            std::string t_contig_name, int t_chrom_index, int t_read_count,
	            std::vector<int> a_five_vec, std::vector<int> b_five_vec,
	            std::vector<int> a_three_vec, std::vector<int> b_three_vec) {

		start = t_start;
		stop = t_stop;
		strand = t_strand;
		chrom_index = t_chrom_index;
		contig_name = t_contig_name;

		five_vec = a_five_vec;
		five_vec.insert(five_vec.end(), b_five_vec.begin(), b_five_vec.end());

		three_vec = a_three_vec;
		three_vec.insert(three_vec.end(), b_three_vec.begin(), b_three_vec.end());

		read_count = t_read_count;
	}

	// Destroy
	~ClusterNode() {
		next = NULL;
		prev = NULL;
	}

	/////////////////////////////////////////////////////////////
	// Get Details
	void set_chrom_index(int t_chrom_index) { chrom_index = t_chrom_index; }
	int get_chrom_index() { return chrom_index; }
	std::string get_contig_name() { return contig_name; }

	int get_strand() { return strand; }

	int get_start() { return start; }
	int get_stop() { return stop; }

	std::string get_headID() { return headID; }
	size_t get_read_count() { return read_count; }

	std::vector<int> get_five_vec() { return five_vec; }
	std::vector<int> get_three_vec() { return three_vec; }
	std::vector<int>* get_five_ref() { return &five_vec; }
	std::vector<int>* get_three_ref() { return &three_vec; }

	std::vector<std::vector<int>>* get_transcripts() { return &transcript_vec; }
	int get_transcript_num() { return transcript_num; }

	/////////////////////////////////////////////////////////////
	// Linking functions
	void set_next(ClusterNode *node) { next = node; }
	void set_prev(ClusterNode *node) { prev = node; }
	ClusterNode* get_next() { return next; }
	ClusterNode* get_prev() { return prev; }

	/////////////////////////////////////////////////////////////
	// Add alignment to cluster
	void add_alignment(int t_5end, int t_3end) {

		five_vec.at(read_count) = t_5end;
		three_vec.at(read_count) = t_3end;
		read_count += 1;

		if (read_count % 1000 == 0) {
			five_vec.resize(five_vec.size() + 1000, 0);
			three_vec.resize(three_vec.size() + 1000, 0);
		}
	}

	// Remove ends of vectors
	void shrink_vectors() {
		five_vec.resize(read_count);
		five_vec.shrink_to_fit();
		three_vec.resize(read_count);
		three_vec.shrink_to_fit();
	}

	// add transcript
	void add_transcript(const std::vector<int> &t_transcript, const float &t_expression) {
		transcript_vec.push_back(t_transcript);
		transcript_expression.push_back(t_expression);
		total_core_points += t_expression;
		transcript_num += 1;
	}

	// Print Transcripts
	void write_transcripts(std::ofstream &gtfFile) {

		int start, stop, x_start, x_stop;
		int regions = 0;
		float prop, quant;
		std::string gene_id = "Unassigned";

		char strand = '+';
		if (this -> get_strand() == 1) { strand = '-'; }

		// gtfFile << "#" << contig_name << ":" << this -> start << "-" << this -> stop
		// 			   << "\t" << read_count << "\t" << this -> five_vec.size() << "\n";

		for (int i = 0; i < transcript_vec.size(); i++) {

			regions	= transcript_vec.at(i).size();

			start = transcript_vec.at(i).at(0);
			stop = transcript_vec.at(i).at(regions - 1);

			// Get expression
			prop = transcript_expression.at(i) / total_core_points;
			if (prop == 1) {
				quant = read_count;
			} else {
				quant = std::round((prop * read_count) * 1000.0f) / 1000.0f;
			}

			// Print Transcript Line
			gtfFile << contig_name << "\timpaqt\ttranscript\t"
			        << start << "\t" << stop << "\t.\t"
			        << strand << "\t.\t"
			        << "gene_id \"" << gene_id << "\";"
			        << " transcript_id \"impaqt."
			        << contig_name << ":"
			        << start << "-" << stop << "\";"
			        << " exons \"" << (regions / 2) << "\";"
			        << " counts \"" << quant << "\";\n";

			// Print Exon Line
			for (int j = 0; j < regions; j += 2) {

				x_start = transcript_vec.at(i).at(j);
				x_stop = transcript_vec.at(i).at(j + 1);

				gtfFile << contig_name << "\timpaqt\texon\t"
				        << x_start << "\t" << x_stop << "\t.\t"
				        << strand  << "\t.\t"
				        << "gene_id \"" << gene_id << "\";"
				        << " transcript_id \"impaqt."
				        << contig_name << ":"
				        << start << "-" << stop << "\";"
				        << " exon \"" << (j / 2) << "\";\n";
			}
		}
	}
};

