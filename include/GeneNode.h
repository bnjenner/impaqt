#pragma once

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Gene Node Class (node in a doubly linked list) */

class GeneNode {

private:

	// Node Details
	std::string geneID;                          // read ID of first read in cluster
	std::string chrom;                           // chromosome index
	int strand = -1;                             // standedness
	int start;                                   // beginning of window
	int stop;                                    // end of window
	int exons = 0;                               // number of exons (or features)
	long double read_count = 0;                  // number of associated reads
	std::vector<int> exon_vec = {0, 0};          // vector for bounds

	// Links
	GeneNode *next = NULL;                       // next ClusterNode
	GeneNode *prev = NULL;                       // pevsious ClusterNode

	/////////////////////////////////////////////////////////////
	/* Private Gene Functions */

	// If two exons overlap
	bool overlap(const int &e1, const int &e2, const int &t1, const int &t2) {
		if (e1 > t2 || e2 < t1) { return false; }
		return true;
	}

	// Insert exon
	void insert_exon(const int &t1, const int &t2) {
		int new_size = exon_vec.size() + 2;
		exon_vec.resize(new_size);
		exon_vec.at(new_size - 2) = t1;
		exon_vec.at(new_size - 1) = t2;
		std::sort(exon_vec.begin(), exon_vec.end()); // sort in ascending order
		exons += 1;
	}

	// Join exons
	void close_gap(const int &t_start, const int &t_stop) {

		// Just populate new vector with exons that are unique
		std::vector<int> t_vec;
		t_vec.reserve(exons * 2); // New exon vec

		for (int i = 0; i < exons; i++) {

			if (exon_vec[(2*i) + 1] < t_start) {
				t_vec.emplace_back(exon_vec[(2*i)]);
				t_vec.emplace_back(exon_vec[(2*i) + 1]);
			
			} else if (exon_vec[(2*i)] > t_stop) {
				t_vec.emplace_back(exon_vec[(2*i)]);
				t_vec.emplace_back(exon_vec[(2*i) + 1]);

			} else {
				t_vec.emplace_back(std::min(exon_vec[(2*i)], t_start));
				t_vec.emplace_back(std::max(exon_vec[(2*i) + 1], t_stop));

				// close gap
				int counter = 0;
				for (int j = i + 1; j < exons; j++) {
					if (exon_vec[(2*j)] < t_vec[t_vec.size() - 1]) {
						t_vec[t_vec.size() - 1] = std::max(t_vec[t_vec.size() - 1], exon_vec[(2*j) + 1]);
						++counter;
					} else { break; }
				}
				i += counter;
			}
		}

		t_vec.shrink_to_fit(); // shrink to fit 
		exon_vec = t_vec; // update exon vector
		exons = exon_vec.size() / 2;
	}


public:

	/////////////////////////////////////////////////////////////
	/* Constructors */

	// Empty
	GeneNode() {}

	// Initialize
	GeneNode(const std::string &geneID, const std::string &chrom, const std::string &strand,
	         const std::string &start, const std::string &stop) {
		
		// Set info
		this -> geneID = geneID;
		this -> chrom = chrom;
		this -> start = std::stoi(start) - 1;
		this -> stop = std::stoi(stop) - 1;
		
		if (strand == "+") {
			this -> strand = 0;
		} else {
			this -> strand = 1;
		}

		exon_vec[0] = this -> start;
		exon_vec[1] = this -> stop;
		exons += 1;
	}

	// Destroy
	~GeneNode() {
		next = NULL;
		prev = NULL;
	}


	/////////////////////////////////////////////////////////////
	/* Get Functions */

	std::string get_chrom() { return chrom; }
	std::string get_geneID() { return geneID; }
	
	int get_strand() { return strand; }
	int get_start() { return start; }
	int get_stop() { return stop; }
	int get_exon_num() { return exons; }
	
	float get_read_count() { return read_count; }

	std::vector<int> get_exon_vec() { return exon_vec; }
	std::vector<int>* get_exon_ref() { return &exon_vec; }

	GeneNode* get_next() { return next; }
	GeneNode* get_prev() { return prev; }

	/////////////////////////////////////////////////////////////
	/* Set functions */

	void set_next(GeneNode *node) { next = node; }
	void set_prev(GeneNode *node) { prev = node; }

	
	/////////////////////////////////////////////////////////////
	/* Gene Functions */

	// Add exon to exon vector
	void add_region(const std::string &str_start, const std::string &str_stop) {

		const int n = this -> get_exon_num();
		const int t_start = std::stoi(str_start) - 1;
		const int t_stop = std::stoi(str_stop) - 1;

		// For all recorded exons
		for (int i = 0; i < n; i++) {

			// If intersection
			if (overlap(exon_vec[(2 * i)], exon_vec[(2 * i) + 1], t_start, t_stop)) {
				exon_vec[(2 * i)] = std::min(exon_vec[(2 * i)], t_start);
				exon_vec[(2 * i) + 1] = std::max(exon_vec[(2 * i) + 1], t_stop);

				// if not last exon and gap now closed
				if (i != (n - 1) && exon_vec[(2 * i) + 1] >= exon_vec[(2 * i) + 2]) {
					close_gap(exon_vec[(2 * i)], exon_vec[(2 * i) + 1]);  // exon spans previous intron
				}

				start = exon_vec[0]; // update start position
				stop = exon_vec.back(); // update stop position
				return;
			}
		}

		insert_exon(t_start, t_stop);
		start = exon_vec[0]; // update start position
		stop = exon_vec.back(); // update stop position
	}

	// Add expression
	void add_expression(const long double &expr) { read_count += expr; }
};


