//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Cluster Node Class (node in a doubly linked list)
class ClusterNode {

private:

	// Node Details
	std::string contig_name;				// name of contig
	int chrom_index;					// chromosome number in index
	int strand = -1;					// standedness
	int start;						// beginning of window
	int stop;						// end of window

	// Read Details
	std::string headID;					// read ID of first read in cluster
	size_t read_count = 0;					// number of associated reads
	size_t vec_count = 0;
	float total_core_points = 0.0;				// number of total core points
	std::vector<int> five_vec;				// vector for 5' ends
	std::vector<int> three_vec;				// vector for 3' ends

	// Links
	ClusterNode *next = NULL;				// next ClusterNode
	ClusterNode *prev = NULL;				// pevsious ClusterNode

	// Transcript Results
	size_t transcript_num = 0;				// number of transcripts identified
	std::vector<std::vector<int>> transcript_vec;		// vector of transcript regions
	std::vector<float> transcript_expression;		// vector of transcript expression
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
	ClusterNode(ClusterNode curr_node, ClusterNode next_node) {

		std::vector<int> tmp_vec; // For copies

		start = curr_node.get_start();
		stop = next_node.get_stop();
		strand = curr_node.get_strand();
		contig_name = curr_node.get_contig_name();
		chrom_index = curr_node.get_chrom_index();

		five_vec = curr_node.get_five_vec();
		tmp_vec = next_node.get_five_vec();
		five_vec.insert(five_vec.end(), tmp_vec.begin(), tmp_vec.end());

		three_vec = curr_node.get_three_vec(); 
		tmp_vec = next_node.get_three_vec();
		three_vec.insert(three_vec.end(), tmp_vec.begin(), tmp_vec.end());

		read_count = curr_node.get_read_count() + next_node.get_read_count();
		vec_count = curr_node.get_vec_count() + next_node.get_vec_count(); 
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
	size_t get_vec_count() { return vec_count; }

	std::vector<int> get_five_vec() { return five_vec; }
	std::vector<int> get_three_vec() { return three_vec; }
	std::vector<int>* get_five_ref() { return &five_vec; }
	std::vector<int>* get_three_ref() { return &three_vec; }
	
	std::vector<std::vector<int>>* get_transcripts() { return &transcript_vec; }
	float get_transcript_expr(const int i) { return transcript_expression.at(i); }
	size_t get_transcript_num() { return transcript_num; }

	/////////////////////////////////////////////////////////////
	// Linking functions
	void set_next(ClusterNode *node) { next = node; }
	void set_prev(ClusterNode *node) { prev = node; }
	ClusterNode* get_next() { return next; }
	ClusterNode* get_prev() { return prev; }

	/////////////////////////////////////////////////////////////
	// Add alignment to cluster
	void add_alignment(const std::vector<int> &positions) {

		// Resize count vecs if necessary
		int n = positions.size();
		for (int i = 1; i <= n; i++) {
			if ((vec_count + i) % 1000 == 0) { 
				five_vec.resize(five_vec.size() + 1000, 0);
				three_vec.resize(three_vec.size() + 1000, 0);
				break;
			}
		}

		// Add positions to 3 and 5 vec
		for (int i = 0; i < n; i++) {
			if (i % 2 == 0) {
				five_vec.at(vec_count) = positions[i];
			} else {
				three_vec.at(vec_count) = positions[i];
				vec_count += 1;
			}
		}

		read_count += 1;
	}

	// Remove ends of vectors
	void shrink_vectors() {
		five_vec.resize(vec_count);
		five_vec.shrink_to_fit();
		three_vec.resize(vec_count);
		three_vec.shrink_to_fit();
	}

	// Sort Point Vectors By 5' Positions
	void sort_vectors() {
		std::vector<int> indices(vec_count);
		std::iota(indices.begin(), indices.end(), 0);
		std::sort(indices.begin(), indices.end(),
			       [&](int i, int j) -> bool {
			            return five_vec[i] < five_vec[j];
			        });

		std::vector<int> tmp_vec(vec_count);
		for (int i = 0; i < vec_count; i++) { tmp_vec[i] = five_vec[indices[i]]; }
		five_vec = tmp_vec;

		for (int i = 0; i < vec_count; i++) { tmp_vec[i] = three_vec[indices[i]]; }
		three_vec = tmp_vec;
	}



	/////////////////////////////////////////////////////////////
	// add transcript
	void add_transcript(const std::vector<int> &t_transcript, const int &t_expression) {
		transcript_vec.push_back(t_transcript);
		transcript_expression.push_back((float)t_expression);
		transcript_assignments.push_back("Unassigned");
		total_core_points += t_expression;
		transcript_num += 1;
	}

	// determine transcript abundance
	void quantify_transcripts() {
		float prop, quant;
		for (int i = 0; i < transcript_vec.size(); i++) {
			prop = transcript_expression.at(i) / total_core_points;
			if (prop == 1) {
				quant = read_count;
			} else {
				quant = std::round((prop * read_count) * 1000.0f) / 1000.0f;
			}
			transcript_expression.at(i) = quant;
		}
	}

	// assign transcripts
	void assign_transcript(const std::string &t_gene_id, const int &i) { transcript_assignments[i] = t_gene_id; }

	// mark as ambiguous
	void assign_ambiguous(const int &i) { transcript_assignments[i] = "Ambiguous"; }

	// Print Transcripts
	void write_transcripts(std::ofstream &gtfFile) {

		int start, stop, x_start, x_stop;
		int regions = 0;
		float quant;
		std::string gene_id;

		// Set Strand
		char strand = '+';
		if (this -> get_strand() == 1) { strand = '-'; }

		for (int i = 0; i < transcript_vec.size(); i++) {

			gene_id = transcript_assignments.at(i);
			regions	= transcript_vec.at(i).size();
			start = transcript_vec.at(i).at(0);
			stop = transcript_vec.at(i).at(regions - 1);
			quant = transcript_expression.at(i);

			// Print Transcript Line
			gtfFile << contig_name << "\timpaqt\ttranscript\t"
			        << start << "\t" << stop << "\t.\t"
			        << strand << "\t.\t"
			        << "gene_id \"" << gene_id << "\";"
			        << " transcript_id \"impaqt."
			        << contig_name << ":"
			        << start << "-" << stop << "\";"
			        << " region \"" << contig_name << ":" 
			        << this -> get_start() << "-" << this -> get_stop() << "\";"
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
				        << " region \"" << contig_name << ":" 
				        << this -> get_start() << "-" << this -> get_stop() << "\";"
				        << " exon \"" << (j / 2) << "\";\n";
			}
		}
	}
};

