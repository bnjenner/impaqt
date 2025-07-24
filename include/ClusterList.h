#pragma once

#include <sstream>
#include <api/BamAux.h>
#include <api/BamReader.h>

#include "global_args.h"
#include "utils.h"
#include "ClusterNode.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Cluster Class (really just a doubly linked list) */

class ClusterList {

private:

	// List Details
	std::string contig_name;
	int contig_index;
	int contig_length = 0;
	int window_size;

	// Links
	ClusterNode *pos_head = NULL;
	ClusterNode *pos_tail = NULL;
	ClusterNode *neg_head = NULL;
	ClusterNode *neg_tail = NULL;

	// Summary
	long double assigned_reads = 0.0;         // Assigned Transcript counts
	long double ambiguous_reads = 0.0;        // Unassigned Transcript counts
	long double unassigned_reads = 0.0;       // Ambigous Transcript counts
	size_t assigned_singles = 0;              // Assigned Read counts
	size_t ambiguous_singles = 0;             // Unssigned Read counts
	size_t unassigned_singles = 0;            // Ambigous Read counts

	size_t multimapped_reads = 0;             // Multimapped Read counts
	size_t low_quality_reads = 0;             // Low Quality Read counts

	size_t total_reads = 0;
	size_t passing_pos_reads = 0;             // Reads passing read check on +
	size_t passing_neg_reads = 0;             // Reads passing read check on -
	size_t transcript_num = 0;

	/////////////////////////////////////////////////////////////
	/* Private Node Methods */
	//void initialize_strand(ClusterNode *&head, ClusterNode *&tail, const int strand, const int &zones);
	void initialize_list(const int t_strand, const int t_pos);
	ClusterNode* extend_list(ClusterNode *&curr, const int t_strand, const int pos);
	void merge_nodes(ClusterNode *&c_node, ClusterNode *&t_head, ClusterNode *&t_tail);
	void delete_list();

public:

	/////////////////////////////////////////////////////////////
	/* Constructors */
	ClusterList() {}
	ClusterList(const int contig_index, const std::string contig_name, const int contig_length) {
		this -> contig_index = contig_index;
		this -> contig_name = contig_name;
		this -> contig_length = contig_length;
		this -> window_size = ImpaqtArguments::Args.window_size;
	}
	~ClusterList() { this -> delete_list(); }

	/////////////////////////////////////////////////////////////
	/* Get Functions */

	std::string get_contig_name() { return contig_name; }

	// Gets
	ClusterNode* get_head(int t_strand) {
		if (t_strand == 0) { return pos_head; }
		return neg_head;
	}
	ClusterNode* get_tail(int t_strand) {
		if (t_strand == 0) { return pos_tail; }
		return neg_tail;
	}

	// Sets
	void set_head(ClusterNode* node, const int t_strand) {
		if (t_strand == 0) {
			pos_head = node;
		} else {
			neg_head = node;
		}
	}
	void set_tail(ClusterNode* node, const int t_strand) {
		if (t_strand == 0) {
			pos_tail = node;
		} else {
			neg_tail = node;
		}
	}

	// Get First cluster by position (just trust me on this one)
	ClusterNode* get_first_cluster(bool &strand) {
		if (pos_head == NULL && neg_head != NULL) {
			strand = 1; return neg_head;
		} else if (neg_head == NULL && pos_head != NULL) {
			strand = 0; return pos_head;
		} else {
			if (pos_head -> get_start() <= neg_head -> get_start()) {
				strand = 0; return pos_head;
			} else {
				strand = 1; return neg_head;
			}
		}
	}

	// Get Next cluster by position
	ClusterNode* get_next_cluster(ClusterNode *&c_node, ClusterNode *&a_prev, ClusterNode *&b_prev, bool &strand) {
		
		a_prev = c_node -> get_next();
		if (a_prev == NULL) {
			c_node = b_prev; strand = !strand;

			// If oppostie strand exhausted or this strand first, continue
		} else if (b_prev == NULL || a_prev -> get_start() <= b_prev -> get_start()) {
			c_node = a_prev;

			// If positives are after negatives, switch strands
		} else { c_node = b_prev; strand = !strand; }

		return c_node;
	}

	// Get Reads Stats
	long double get_assigned_reads() { return assigned_reads + (long double)assigned_singles; }
	long double get_unassigned_reads() { return unassigned_reads + (long double)unassigned_singles; }
	long double get_ambiguous_reads() { return ambiguous_reads + (long double)ambiguous_singles; }
	size_t get_multimapped_reads() { return multimapped_reads; }
	size_t get_low_quality_reads() { return low_quality_reads; }
	size_t get_total_reads() { return total_reads; }
	size_t get_passing_reads(const int &strand) {
		if (strand == 0) { return passing_pos_reads; }
		return passing_neg_reads;
	}

	// Get Transcript Number
	size_t get_transcript_num() {
		ClusterNode *node = pos_head;
		for (int i = 0; i < 2; i ++) {
			if (i != 0) { node = neg_head; }
			while (node != NULL) {
				if (!(node -> is_skipped())) { transcript_num += node -> get_transcript_num(); }
				node = node -> get_next();
			}
		}
		return transcript_num;
	}

	/////////////////////////////////////////////////////////////
	/* Private Alignment Methods */
	void calculate_splice(BamTools::BamAlignment &alignment, std::vector<int> &positions);
	bool read_check(const BamTools::BamAlignment &alignment);

	/////////////////////////////////////////////////////////////
	/* Counting Functions */
	void add_assigned_reads(const long double &expr) { assigned_reads += expr; }
	void add_unassigned_reads(const long double &expr) { unassigned_reads += expr; }
	void add_ambiguous_reads(const long double &expr) { ambiguous_reads += expr; }
	void add_assigned_singles(const size_t &expr) { assigned_reads += expr; }
	void add_unassigned_singles(const size_t &expr) { unassigned_reads += expr; }
	void add_ambiguous_singles(const size_t &expr) { ambiguous_reads += expr; }

	/////////////////////////////////////////////////////////////
	/* List Functions */

	// Find Nearest Region in List
	void jump_to_cluster(ClusterNode *&node, const int &pos) {
		while (node -> get_next() != NULL) {
			if (!(node -> read_contained(pos))) {
				node = node -> get_next();
			} else { break; }
		}
	}

	/////////////////////////////////////////////////////////////
	/* Cluster Functions */
	bool create_clusters(BamTools::BamReader &inFile, BamTools::BamAlignment &alignment);
	void collapse_clusters(int t_strand);

	/////////////////////////////////////////////////////////////
	/* Output Functions */
	void write_clusters_as_GTF(std::ofstream &gtfFile);

	/////////////////////////////////////////////////////////////
	/* Functions for Testing Suite */
	void print_clusters(int t_strand) {
		ClusterNode *node = get_head(t_strand);
		while (node != NULL) {
			std::cout << node -> get_contig_name() << "\t"
			          << node -> get_start() << "\t" << node -> get_stop() << "\t"
			          << node -> get_read_count() << "\n";
			node = node -> get_next();
		}
	}

	std::string string_clusters(int t_strand) {
		std::stringstream ss;
		ClusterNode *node = get_head(t_strand);
		while (node != NULL) {
			ss << node -> get_contig_name() << "\t"
			   << node -> get_start() << "\t" << node -> get_stop() << "\t"
			   << node -> get_read_count() << "\n";
			node = node -> get_next();
		}
		return ss.str();
	}
};
