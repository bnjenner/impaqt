#pragma once

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* DBSCAN and Related Functions */

// Links a 5' DBSCAN cluster to a 3' DBSCAN cluster. -1 means that prime is
// unassigned. Ordered by {five, three} with -1 sorting first, which matches the
// iteration order of the 2-char-string keys this replaced.
struct Path {
	int five = -1;
	int three = -1;
	bool operator<(const Path &o) const {
		if (five != o.five) { return five < o.five; }
		return three < o.three;
	}
};

// Check if transcripts overlap / are contained in another transcript
bool check_subset(const std::vector<int>& a, const std::vector<int>& b);

// Get Core Points of Transcript
int get_quant(const std::vector <int> result, const std::vector<std::vector<int>> &init_copy, const std::vector<int> counts);

// Report Unique Transcripts (no overlapping, used for Mitochrondria)
void report_transcripts(ClusterNode *node, std::vector<std::vector<int>> &result, std::vector<int> &counts);

// Merge Overlapping Transcripts
bool overlap_aux(std::vector<std::vector<int>> &transcripts);

// Reduce Transcript Number by Overlapping. Report Unique Transcripts
void overlap_clusters(ClusterNode *curr_node, std::vector<std::vector<int>> &transcripts, std::vector<int> &counts);

/// Get Transcript Coordinates
void get_coordinates(const std::map<Path, int> &paths,
                     const std::map<int, std::vector<int>> &regions_5, const std::map<int, std::vector<int>> &regions_3,
                     std::vector<std::vector<int>> *transcripts, std::vector<int> *counts);

// Find all linked DBSCAN clusters
void get_linked_clusters(ClusterNode *curr_node, std::map<Path, int> &path_map,
                         const std::vector<int> &assign_5, const std::vector<int> &assign_3);

// Get Nearest neighbors for Indices in DBSCAN
std::vector<int> get_nearest_neighbors(const int &i, const int &points, std::vector<bool> &queued,
                                       const std::vector<int> &indices, const std::vector<int> *adj_vec);


// DBSCAN Clustering Function, inspired by https://github.com/Eleobert/dbscan/blob/master/dbscan.cpp
std::vector<int> dbscan(ClusterNode *curr_node, const int &points, const int &min_counts,
                        std::map<int, std::vector<int>> &regions, const bool &five);

// Initiate Transcript Identifying Procedure
void identify_transcripts_dbscan(ClusterList *cluster,  const int &strand);

// Merge identified transscripts
void merge_ranscripts(ClusterNode *c_node, ClusterNode *n_node);

// Combine clusters with nonzero neighbors
void collapse_transcripts(ClusterList *cluster, int t_strand);

// Combine clusters with nonzero neighbors
void identify_transcripts(ClusterList *cluster, int t_strand);