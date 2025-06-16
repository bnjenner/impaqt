#include "GeneNode.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Annotation Class (doubley linked list)
class AnnotationList {

private:

	std::string annotation_file;
	std::string stranded;
	bool isGFF;

	// Get First Feature of Chrom
	std::unordered_map<std::string, GeneNode*> pos_chrom_map;
	std::unordered_map<std::string, GeneNode*> neg_chrom_map;

	int features = 0;
	std::string feature_id;
	std::string feature_tag;

	GeneNode *pos_head = NULL;
	GeneNode *pos_tail = NULL;
	GeneNode *neg_head = NULL;
	GeneNode *neg_tail = NULL;

	/////////////////////////////////////////////////////////////
	// Create new gene node
	GeneNode* create_new_node(std::vector<std::string> &columns) {
		GeneNode *temp = new GeneNode(columns[8], 	// Feature ID
		                              columns[0], 	// Chrom
		                              columns[6],	// Strand
		                              columns[3],	// Start
		                              columns[4]);	// Stop
		return temp;
	}

	void set_pos_head(std::vector<std::string> &columns) {
		pos_head = create_new_node(columns);
		pos_tail = pos_head;
	}

	void set_neg_head(std::vector<std::string> &columns) {
		neg_head = create_new_node(columns);
		neg_tail = neg_head;
	}

	void pos_extend(std::vector<std::string> &columns) {
		// Create new gene node
		GeneNode *temp = create_new_node(columns);
		pos_tail -> set_next(temp);
		temp -> set_prev(pos_tail);
		pos_tail = temp;
	}

	void neg_extend(std::vector<std::string> &columns) {
		GeneNode *temp = create_new_node(columns);
		neg_tail -> set_next(temp);
		temp -> set_prev(neg_tail);
		neg_tail = temp;
	}

	/////////////////////////////////////////////////////////////
	// Get feature ID
	bool set_feature_id(std::vector<std::string> &t_columns) {

		int i;
		char sep = '"';
		std::string tag;
		std::string value;

		// if gff, different separator.
		if (ImpaqtArguments::Args.isGFF) { sep = '='; }

		// Iterate through annotation features
		std::istringstream anno_stream(t_columns.at(8));
		while (std::getline(anno_stream, tag, ';')) {

			if (tag.find(feature_id) == std::string::npos) { continue; }

			i = 0;
			std::istringstream feature_stream(tag);
			while (std::getline(feature_stream, value, sep)) {
				if (i % 2 == 1) { t_columns[8] = value; return true; }
				i++;
			}
		}

		return false;
	}


public:

	/////////////////////////////////////////////////////////////
	// Proper Constructor
	AnnotationList() {
		annotation_file = ImpaqtArguments::Args.annotation_file;
		feature_id = ImpaqtArguments::Args.feature_id;
		feature_tag = ImpaqtArguments::Args.feature_tag;
		stranded = ImpaqtArguments::Args.stranded;
		isGFF = ImpaqtArguments::Args.isGFF;
	}

	// Proper Destructor
	~AnnotationList() {

		GeneNode *curr_node = pos_head;
		GeneNode *temp_node = NULL;

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

	/////////////////////////////////////////////////////////////
	// Get Head Node
	GeneNode* get_head(const int t_strand) {
		if (t_strand == 0) { return pos_head; }
		return neg_head;
	}

	// Get Tail Node
	GeneNode* get_tail(const int t_strand) {
		if (t_strand == 0) { return pos_tail; }
		return neg_tail;
	}

	GeneNode* jump_to_chrom(const std::string t_chrom, const int t_strand) {
		if (t_strand == 0) { 
			if (pos_chrom_map.find(t_chrom) == pos_chrom_map.end()) { return NULL; }
			return pos_chrom_map[t_chrom];
		} else {
			if (neg_chrom_map.find(t_chrom) == neg_chrom_map.end()) { return NULL; }
			return neg_chrom_map[t_chrom];
		}
	}

	/////////////////////////////////////////////////////////////
	// Get Features
	int get_features() { return features; }

	/////////////////////////////////////////////////////////////
	// Create Gene Graph Structure
	void create_gene_list() {

		// Open file
		std::ifstream infile(annotation_file);
		if (!infile) { throw "ERROR: Could not read annotation file."; }

		// Iterate through lines in file
		std::string col;
		std::string line;
		std::vector<std::string> columns{9, ""};
		while (std::getline(infile, line)) {

			// skip headers
			if (line[0] == '#') { continue; }

			// populate column vector
			int i = 0;
			std::istringstream column_stream(line);
			while (std::getline(column_stream, col, '\t')) {
				columns.at(i) = col;
				i += 1;
			}

			// if not feature tag
			if (columns[2] != ImpaqtArguments::Args.feature_tag) { continue; }

			// if feature id not found
			if (!set_feature_id(columns)) {
				std::cerr << "ERROR: Could not find feature tag in line:\n" << line << "\n";
				throw "ERROR: Could not find feature tag in line of annotation file. Check consistency of formatting.";
			}

			features += 1;

			// If positive strand
			if (columns[6] == "+") {

				// If first gene on positive strand
				if (pos_head == NULL) {
					set_pos_head(columns);
				} else {
					// If same gene ID
					if (columns[8] == pos_tail -> get_geneID()) {
						pos_tail -> add_region(columns[3], columns[4]);
					} else {
						pos_extend(columns); // create new gene node
					}
				}	

				// If New Chrom add to chrom map
				if (pos_chrom_map.find(columns[0]) == pos_chrom_map.end()) {
					pos_chrom_map[columns[0]] = pos_tail;
				}

				// If negative strand
			} else if (columns[6] == "-") {

				// If first gene on positive strand
				if (neg_head == NULL) {
					set_neg_head(columns);
				} else {
					// If same gene ID
					if (columns[8] == neg_tail -> get_geneID()) {
						neg_tail -> add_region(columns[3], columns[4]);
					} else {
						neg_extend(columns); // create new gene node
					}
				}

				// If New Chrom add to chrom map
				if (neg_chrom_map.find(columns[0]) == neg_chrom_map.end()) {
					neg_chrom_map[columns[0]] = neg_tail;
				}

			} else {
				std::cerr << "ERROR: Unrecognized symbol in strand column:\n" << line << "\n";
				throw "ERROR: Unrecognized symbol in strand column of annotation file. Check consistency of formatting.";
			}
		}
	}

	/////////////////////////////////////////////////////////////
	// Print genes and counts
	void print_gene_counts() {

		// Return if empty GTF?
		if (pos_head == NULL && neg_head == NULL) {
			std::cerr << "// NOTICE: No genes found in annotation file.\n";
			return;
		}


		bool strand;
		GeneNode *curr_node;
		GeneNode *prev_pos_node = pos_head;
		GeneNode *prev_neg_node = neg_head;

		// Get First Gene
		if (pos_head == NULL && neg_head != NULL) {
			curr_node = neg_head; strand = 1;
		} else if (neg_head == NULL && pos_head != NULL) {
			curr_node = pos_head; strand = 0;
		} else {
			if (pos_head -> get_start() < neg_head -> get_start()) {
				curr_node = pos_head; strand = 0;
			} else {
				curr_node = neg_head; strand = 1;
			}
		}

		// Iterate Throught Genes
		while (true) {

			// Break When all genes have been exhausted
			if (prev_neg_node == NULL && prev_pos_node == NULL) { break; }

			// Report Gene and Counts
			std::cout << curr_node -> get_geneID() << "\t" 
					  << curr_node -> get_read_count() << "\n";
			
			// Strand switching conditions :(
			if (strand == 0) {
				prev_pos_node = curr_node -> get_next();

				// If positives exhausted, switch strands
				if (prev_pos_node == NULL) {
					curr_node = prev_neg_node; strand = 1;
					continue;
				}

				// If negatives exhausted, continue with positives
				if (prev_neg_node == NULL || prev_pos_node -> get_start() < prev_neg_node -> get_start()) {
					curr_node = prev_pos_node; strand = 0;

					// If positives are after negatives, switch strands
				} else {
					curr_node = prev_neg_node; strand = 1;
				}

			} else {
				prev_neg_node = curr_node -> get_next();

				// If negatives exhausted, switch strands
				if (prev_neg_node == NULL) { 
					curr_node = prev_pos_node; strand = 0;
					continue;
				}

				// If positives exhausted, continue with negatives
				if (prev_pos_node == NULL || prev_neg_node -> get_start() < prev_pos_node -> get_start()) {
					curr_node = prev_neg_node; strand = 1;

					// If negatives are after positives, switch strands
				} else {
					curr_node = prev_pos_node; strand = 0;
				}
			}
		}
	}


	/////////////////////////////////////////////////////////////
	// Mainly for Test Suite
	// Print genes into strings for tests
	std::string string_genes(const int &strand) {
		std::stringstream ss;
		GeneNode *curr_node = get_head(strand);
		while (curr_node != NULL) {
			ss << curr_node -> get_chrom() << "\t" << curr_node -> get_geneID() << "\t";
			for (const auto &region : curr_node -> get_exon_vec()) { ss << region << "\t"; }
			ss << "\n";
			curr_node = curr_node -> get_next();
		}
		return ss.str();
	}


	// Print genes 
	void print_genes(const int &strand) {
		GeneNode *curr_node = get_head(strand);
		while (curr_node != NULL) {
			std::cerr << curr_node -> get_chrom() << ": " << curr_node -> get_geneID() << "\n";
			for (int i = 0; i < curr_node -> get_exon_num(); i++) {
				std::cerr << "\t" << curr_node -> get_exon_vec()[(2*i)] << "\t"
								  << curr_node -> get_exon_vec()[(2*i)+1] << "\n";
			}
			curr_node = curr_node -> get_next();
		}
	}

	// Print genes and counts
	void print_chrom_map() {
		std::cerr << "// Positive Strand Chrom Map\n";
		for (const auto &pair : pos_chrom_map) {
			std::cerr << pair.first << ": " << pair.second << "\n";
		}

		std::cerr << "// Negative Strand Chrom Map\n";
		for (const auto &pair : neg_chrom_map) {
			std::cerr << pair.first << ": " << pair.second << "\n";
		}
	
	}
};