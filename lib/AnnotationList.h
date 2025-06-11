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
				if (pos_head == NULL) { set_pos_head(columns); continue; }

				// If same gene ID
				if (columns[8] == pos_tail -> get_geneID()) {
					pos_tail -> add_region(columns[3], columns[4]);
				} else {
					pos_extend(columns); // create new gene node
					
					// If New Chrom add to chrom map
					if (pos_chrom_map.find(columns[0]) == pos_chrom_map.end()) { pos_chrom_map[columns[0]] = pos_tail; }
				}

				// If negative strand
			} else if (columns[6] == "-") {

				// If first gene on positive strand
				if (neg_head == NULL) {  set_neg_head(columns); continue; }

				// If same gene ID
				if (columns[8] == neg_tail -> get_geneID()) {
					neg_tail -> add_region(columns[3], columns[4]);
				} else {
					neg_extend(columns); // create new gene node

					// If New Chrom add to chrom map
					if (neg_chrom_map.find(columns[0]) == neg_chrom_map.end()) { neg_chrom_map[columns[0]] = neg_tail; }
				}

			} else {
				std::cerr << "ERROR: Unrecognized symbol in strand column:\n" << line << "\n";
				throw "ERROR: Unrecognized symbol in strand column of annotation file. Check consistency of formatting.";
			}
		}
	}

	/////////////////////////////////////////////////////////////
	// Print genes and counts
	void print_genes() {
		GeneNode *curr_node = pos_head;
		while (curr_node != NULL) {
			std::cout << curr_node -> get_chrom() << "\t" << curr_node -> get_geneID() << "\t";
			for (const auto &region : curr_node -> get_exon_vec()) { std::cout << region << "\t"; }
			std::cout << "\n";
			curr_node = curr_node -> get_next();
		}
	}

	// Print genes into strings for tests
	std::string string_genes() {
		std::stringstream ss;
		GeneNode *curr_node = pos_head;
		while (curr_node != NULL) {
			ss << curr_node -> get_chrom() << "\t" << curr_node -> get_geneID() << "\t";
			for (const auto &region : curr_node -> get_exon_vec()) { ss << region << "\t"; }
			ss << "\n";
			curr_node = curr_node -> get_next();
		}
		return ss.str();
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