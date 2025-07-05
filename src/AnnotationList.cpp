#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <unordered_map>

#include "global_args.h"
#include "AnnotationList.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Annotation Class (doubley linked list) */

// Get feature ID
bool AnnotationList::set_feature_id(std::vector<std::string> &t_columns) {

	int i; char sep = '"';
	std::string tag, value;

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

void AnnotationList::set_head(const std::vector<std::string> &columns) {
	if (columns[6] == "+") {
		AnnotationList::pos_head = create_new_node(columns);
		AnnotationList::pos_tail = AnnotationList::pos_head;
	} else {
		AnnotationList::neg_head = create_new_node(columns);
		AnnotationList::neg_tail = AnnotationList::neg_head;
	}
}

void AnnotationList::extend(const std::vector<std::string> &columns) {
	GeneNode *gene = create_new_node(columns);
	if (columns[6] == "+") {
		AnnotationList::pos_tail -> set_next(gene);
		gene -> set_prev(AnnotationList::pos_tail);
		AnnotationList::pos_tail = gene;
	} else {
		AnnotationList::neg_tail -> set_next(gene);
		gene -> set_prev(AnnotationList::neg_tail);
		AnnotationList::neg_tail = gene;
	}
}

void AnnotationList::add_line(const std::vector<std::string> &columns) {

	// Get proper strand
	GeneNode **head = &(AnnotationList::pos_head);
	GeneNode **tail = &(AnnotationList::pos_tail);
	std::unordered_map<std::string, GeneNode*> *chrom_map = &(AnnotationList::pos_chrom_map);
	
	if (columns[6] == "-") {
		head = &(AnnotationList::neg_head); tail = &(AnnotationList::neg_tail);
		chrom_map = &(AnnotationList::neg_chrom_map);
	}

	// If first gene on strand
	if ((*head) == NULL) {
		AnnotationList::set_head(columns);

	} else {
		// If same gene ID
		if (columns[8] == (*tail) -> get_geneID()) {
			(*tail) -> GeneNode::add_region(columns[3], columns[4]);
		
		} else { AnnotationList::extend(columns); } // create new gene node
	}	

	// If new chrom add to chrom map
	if (chrom_map -> find(columns[0]) == chrom_map -> end()) {
		(*chrom_map)[columns[0]] = *tail;
	}
}

/////////////////////////////////////////////////////////////
/* List Functions */

// Create Gene Graph Structure
void AnnotationList::create_gene_list() {

	// Open file
	std::ifstream infile(AnnotationList::annotation_file);
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
		if (!AnnotationList::set_feature_id(columns)) {
			std::cerr << "ERROR: Could not find feature tag in line:\n" << line << "\n";
			throw "ERROR: Could not find feature tag in line of annotation file. Check consistency of formatting.";
		}

		AnnotationList::add_line(columns);
		++AnnotationList::features;
	}
}

/////////////////////////////////////////////////////////////
/* Output Functions */

// Print genes and counts
void  AnnotationList::print_gene_counts() {

	// Return if empty GTF?
	if (AnnotationList::pos_head == NULL &&  AnnotationList::neg_head == NULL) {
		std::cerr << "// NOTICE: No genes found in annotation file.\n";
		return;
	}

	bool strand;
	GeneNode *prev_pos = AnnotationList::pos_head;
	GeneNode *prev_neg = AnnotationList::neg_head;
	GeneNode *node = AnnotationList::get_first_gene(strand);

	// Iterate Throught Genes
	while (true) {

		// Break When all genes have been exhausted
		if (prev_neg == NULL && prev_pos == NULL) { break; }

		// Report Gene and Counts
		std::cout << node -> GeneNode::get_geneID() << "\t" 
				  << node -> GeneNode::get_read_count() << "\n";
		
		// Get Next Cluster by position
		if (strand == 0) {
			node = AnnotationList::get_next_gene(node, prev_pos, prev_neg, strand);
		} else {
			node = AnnotationList::get_next_gene(node, prev_neg, prev_pos, strand);
		}
	}
}