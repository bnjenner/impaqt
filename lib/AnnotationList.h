#include "GeneNode.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Annotation Class (doubley linked list) */

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
	/* Private Gene Methods */

	// Create new gene node
	GeneNode* create_new_node(const std::vector<std::string> &columns) {
		GeneNode *temp = new GeneNode(columns[8], 	// Feature ID
		                              columns[0], 	// Chrom
		                              columns[6],	// Strand
		                              columns[3],	// Start
		                              columns[4]);	// Stop
		return temp;
	}

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


	/////////////////////////////////////////////////////////////
	/* Private List Methods */

	void set_head(const std::vector<std::string> &columns) {
		if (columns[6] == "+") {
			pos_head = create_new_node(columns);
			pos_tail = pos_head;
		} else {
			neg_head = create_new_node(columns);
			neg_tail = neg_head;
		}
	}

	void extend(const std::vector<std::string> &columns) {
		GeneNode *t_gene = create_new_node(columns);
		if (columns[6] == "+") {
			pos_tail -> set_next(t_gene);
			t_gene -> set_prev(pos_tail);
			pos_tail = t_gene;
		} else {
			neg_tail -> set_next(t_gene);
			t_gene -> set_prev(neg_tail);
			neg_tail = t_gene;
		}
	}

	void add_line(const std::vector<std::string> &columns) {

		// Get proper strand
		GeneNode **head = &pos_head;
		GeneNode **tail = &pos_tail;
		std::unordered_map<std::string, GeneNode*> *chrom_map = &pos_chrom_map;
		if (columns[6] == "-") {
			head = &neg_head; tail = &neg_tail;
			chrom_map = &neg_chrom_map;
		}

		// If first gene on strand
		if ((*head) == NULL) {
			set_head(columns);
		} else {
			// If same gene ID
			if (columns[8] == (*tail) -> get_geneID()) {
				(*tail) -> add_region(columns[3], columns[4]);
			} else {
				extend(columns); // create new gene node
			}
		}	

		// If new chrom add to chrom map
		if (chrom_map -> find(columns[0]) == chrom_map -> end()) {
			(*chrom_map)[columns[0]] = *tail;
		}

	}

public:

	/////////////////////////////////////////////////////////////
	/* Constructors */

	// Proper
	AnnotationList() {
		annotation_file = ImpaqtArguments::Args.annotation_file;
		feature_id = ImpaqtArguments::Args.feature_id;
		feature_tag = ImpaqtArguments::Args.feature_tag;
		stranded = ImpaqtArguments::Args.stranded;
		isGFF = ImpaqtArguments::Args.isGFF;
	}

	// Destrpy
	~AnnotationList() {
		GeneNode *c_node = pos_head;
		GeneNode *t_node = NULL;
		for (int i = 0; i < 2; i ++) {
			if (i != 0) { c_node = neg_head; t_node = NULL; }
			while (c_node != NULL) {
				t_node = c_node;
				c_node = c_node -> get_next();
				delete t_node;
			}
		}
	}


	/////////////////////////////////////////////////////////////
	/* Get Functions */

	int get_features() { return features; }

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

	// Get first gene by position (just trust me on this one)
	GeneNode* get_first_gene(GeneNode *pos, GeneNode *neg, bool &strand) {
		if (pos == NULL && neg != NULL) {
			strand = 1; return neg;
		} else if (neg == NULL && pos != NULL) {
			strand = 0; return pos;
		} else {
			if (pos -> get_start() < neg -> get_start()) {
				strand = 0; return pos;
			} else {
				strand = 1; return neg;
			}
		}
	}

	// Get next gene by position
	GeneNode* get_next_gene(GeneNode *&c_node, GeneNode *&a_prev, GeneNode *&b_prev, bool &strand) {
		
		a_prev = c_node -> get_next();
		if (a_prev == NULL) {
			c_node = b_prev; strand = !strand;

			// If oppostie strand exhausted or this strand first, continue
		} else if (b_prev == NULL || a_prev -> get_start() < b_prev -> get_start()) {
			c_node = a_prev;

			// If positives are after negatives, switch strands
		} else { c_node = b_prev; strand = !strand; }

		return c_node;
	}


	/////////////////////////////////////////////////////////////
	/* List Functions */

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

			add_line(columns);
			features += 1;
		}
	}


	/////////////////////////////////////////////////////////////
	/* Output Functions */

	// Print genes and counts
	void print_gene_counts() {

		// Return if empty GTF?
		if (pos_head == NULL && neg_head == NULL) {
			std::cerr << "// NOTICE: No genes found in annotation file.\n";
			return;
		}

		bool strand;
		GeneNode *prev_pos_node = pos_head;
		GeneNode *prev_neg_node = neg_head;
		GeneNode *c_node = get_first_gene(pos_head, neg_head, strand) ;

		// Iterate Throught Genes
		while (true) {

			// Break When all genes have been exhausted
			if (prev_neg_node == NULL && prev_pos_node == NULL) { break; }

			// Report Gene and Counts
			std::cout << c_node -> get_geneID() << "\t" 
					  << c_node -> get_read_count() << "\n";
			
			// Get Next Cluster by position
			if (strand == 0) {
				c_node = get_next_gene(c_node, prev_pos_node, prev_neg_node, strand);
			} else {
				c_node = get_next_gene(c_node, prev_neg_node, prev_pos_node, strand);
			}
		}
	}


	/////////////////////////////////////////////////////////////
	/* Cluster Functions */

	// Print genes into strings for tests
	std::string string_genes(const int &strand) {
		std::stringstream ss;
		GeneNode *c_node = get_head(strand);
		while (c_node != NULL) {
			ss << c_node -> get_chrom() << "\t" << c_node -> get_geneID() << "\t";
			for (const auto &region : c_node -> get_exon_vec()) { ss << region << "\t"; }
			ss << "\n";
			c_node = c_node -> get_next();
		}
		return ss.str();
	}


	// Print genes 
	void print_genes(const int &strand) {
		GeneNode *c_node = get_head(strand);
		while (c_node != NULL) {
			std::cerr << c_node -> get_chrom() << ": " << c_node -> get_geneID() << "\n";
			for (int i = 0; i < c_node -> get_exon_num(); i++) {
				std::cerr << "\t" << c_node -> get_exon_vec()[(2*i)] << "\t"
								  << c_node -> get_exon_vec()[(2*i)+1] << "\n";
			}
			c_node = c_node -> get_next();
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