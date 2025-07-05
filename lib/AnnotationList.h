#pragma once

#include <sstream>
#include <algorithm>
#include <unordered_map>

#include "global_args.h"
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
	bool set_feature_id(std::vector<std::string> &t_columns);

	/////////////////////////////////////////////////////////////
	/* Private List Methods */

	void set_head(const std::vector<std::string> &columns);
	void extend(const std::vector<std::string> &columns);
	void add_line(const std::vector<std::string> &columns);

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
	GeneNode* get_first_gene(bool &strand) {
		if (pos_head == NULL && neg_head != NULL) {
			strand = 1; return neg_head;
		} else if (neg_head == NULL && pos_head != NULL) {
			strand = 0; return pos_head;
		} else {
			if (pos_head -> get_start() < neg_head -> get_start()) {
				strand = 0; return pos_head;
			} else {
				strand = 1; return neg_head;
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
	void create_gene_list();

	/////////////////////////////////////////////////////////////
	/* Output Functions */

	// Print genes and counts
	void print_gene_counts();

	/////////////////////////////////////////////////////////////
	/* Functions For Test Suite */

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