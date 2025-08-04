#include <iostream>

#include "ClusterList.h"
#include "DBSCAN.h"
#include "ContainmentList.h"
#include "utils.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* DBSCAN and Related Functions */


// Get Core Points of Transcript
int get_quant(const std::vector <int> result, const std::vector<std::vector<int>> &init_copy, const std::vector<int> counts) {
	int core_points = 0;
	for (int i = 0; i < init_copy.size(); i++) {
		if (check_containment_strict(init_copy[i], result)) {
			core_points += counts[i];
		}
	}
	return core_points;
}

// Report Unique Transcripts (no overlapping, used for Mitochrondria)
void report_transcripts(ClusterNode *node, std::vector<std::vector<int>> &result, std::vector<int> &counts) {
	for (int i = 0; i < result.size(); i++) { node -> add_transcript(result[i], counts[i]); }
}


// Merge Overlapping Transcripts
bool overlap_aux(std::vector<std::vector<int>> &transcripts) {

	// Construct containment list 
	ContainmentList head;
	ContainmentList *curr;
	std::vector<int> visited = {};

	// Construct containment list
	bool unique = true;
	for (int i = 0; i < transcripts.size(); i++) {

		// Set head
		if (head.indices == 0) {
			head = ContainmentList(transcripts[i]);
			curr = &head;

			// or extend list
		} else if (std::find(visited.begin(), visited.end(), i) == visited.end()) {
			curr -> set_next(new ContainmentList(transcripts[i]));
			curr = curr -> next;

		} else { continue; }

		for (int j = i + 1; j < transcripts.size(); j++) {

			// Skip if overlap notpossible
			if (transcripts[j][0] > curr -> get_back()) { continue; }

			if (check_containment(transcripts[j], transcripts[i])) {
				curr -> add_interval(transcripts[j]);
				visited.push_back(j);
				unique = false;
			}
		}
		curr -> collapse_intervals();
	}


	curr = &head;
	transcripts.clear();
	while (curr != nullptr) {
		// Add if unique
		if (curr -> sublist_count == 0) { 
			transcripts.push_back(curr -> vals); 

			// Add if extended
		} else {
			for (const auto &n : curr -> sublist) { transcripts.push_back(n -> vals); }
		}

		curr -> clean();
		curr = curr -> next;
	}

	return unique;
}


// Reduce Transcript Number by Overlapping. Report Unique Transcripts
void overlap_clusters(ClusterNode *curr_node, std::vector<std::vector<int>> &transcripts, std::vector<int> &counts) {

	std::vector<std::vector<int>> init_copy = transcripts;

	// Return if no need to overlap
	if (transcripts.size() == 1) { return; }

	// Reverse and Negative if reverse strand
	if (curr_node -> get_strand() == 1) {
		reverse_transcripts(transcripts);
	} else {
		std::sort(transcripts.begin(), transcripts.end(), compare_first_element);
	}

	bool unique = overlap_aux(transcripts);
	while (!unique) {
		unique = overlap_aux(transcripts);
	}

	// Reverse and Negative Results if Necessary
	if (curr_node -> get_strand() == 1) {
		reverse_transcripts(transcripts);
	} else {
		std::sort(transcripts.begin(), transcripts.end(), compare_first_element);
	}
	
	// If more than one transcript identified
	int core_points = 0;
	if (transcripts.size() == 1) {
		for (const auto &c : counts) { core_points += c; }
		counts = {core_points};
	} else {
		std::vector<int> new_counts(transcripts.size(), 0);
		for (int i = 0; i < transcripts.size(); i++) {
			core_points = get_quant(transcripts[i], init_copy, counts);
			new_counts[i] = core_points;
		}
		counts = new_counts;
	}
}


// Get Transcript Coordinates
void get_coordinates(ClusterNode *curr_node, const std::map<std::string, int> &paths,
                     std::vector<std::vector<int>> &core_5, std::vector<std::vector<int>> &core_3,
                     std::vector<std::vector<int>> *transcripts, std::vector<int> *counts) {

	// Get transcript coordinates
	//	BNJ: 5/2/2025 - Horrible extending vector, but they're short so it shouldn't matter THAT much
	//	BNJ: 5/31/2025 - Also worth mentioning, the tmp_vec should never be more than 4 in length
	//	BNJ: 6/4/2025 - Man, I should really use bits instead of a string here.

	int index, min_pos, max_pos; 

	for (const auto &p : paths) {

		if (p.second == 0) { continue; }

		std::vector<int> tmp_vec;

		// add 5' region
		if (p.first[0] != '-') {
			index = std::stoi(p.first.substr(0, 1));
			min_pos = get_pos_min(index, core_5, curr_node -> get_five_ref());
			max_pos = get_pos_max(index, core_5, curr_node -> get_five_ref());
			if (min_pos > max_pos) { variable_swap(min_pos, max_pos); } 
			tmp_vec.push_back(min_pos);
			tmp_vec.push_back(max_pos);
		}

		// add 3' region
		if (p.first[1] != '-') {
			index = std::stoi(p.first.substr(1, 1));
			min_pos = get_pos_min(index, core_3, curr_node -> get_three_ref());
			max_pos = get_pos_max(index, core_3, curr_node -> get_three_ref());
			if (min_pos > max_pos) { variable_swap(min_pos, max_pos); } 
			tmp_vec.push_back(min_pos);
			tmp_vec.push_back(max_pos);
		}

		// If two regions and they are close or out of order, merge
		if (tmp_vec.size() > 2) {
			if (ImpaqtArguments::Args.epsilon >= std::abs(tmp_vec[2] - tmp_vec[1])) {
				tmp_vec = {tmp_vec[0], tmp_vec[3]};
			} else if (tmp_vec[2] <= tmp_vec[1]) {
				tmp_vec = {std::min(tmp_vec[0], tmp_vec[2]), std::max(tmp_vec[1], tmp_vec[3])};
			}
		}

		transcripts -> emplace_back(std::move(tmp_vec));
		counts -> push_back(p.second);
	}
}


// Find all linked DBSCAN clusters
void get_linked_clusters(ClusterNode *curr_node, std::map<std::string, int> &path_map,
                         const std::vector<int> &assign_5, const std::vector<int> &assign_3) {

	std::string path;

	for (int i = 0; i < curr_node -> get_vec_count(); i++) {

		path = "";

		// assigned in 5' DBSCAN
		if (assign_5.at(i) != -1) {

			path = std::to_string(assign_5.at(i)) + '-';
			if (assign_3.at(i) != -1) { path.at(1) = std::to_string(assign_3.at(i))[0]; }

			// unassigned in 5' DBSCAN
		} else if (assign_5.at(i) == -1 && assign_3.at(i) != -1) {
			path = "-" + std::to_string(assign_3.at(i));
		}

		if (path == "") { continue; }

		if (path_map.find(path) != path_map.end()) {
			path_map[path] += 1; // Increment count if path already exists
		} else {
			path_map[path] = 1; // Add new path with count of 1
		}
	}


	// Absorb orphan paths
	int pos;
	for (const auto& p1 : path_map) {

		if (p1.second == 0) { continue; }

		// if path is orphaned
		pos = p1.first.find('-');
		if (pos != std::string::npos) {
			
			for (const auto& p2 : path_map) { 

				if (p1.first == p2.first) { continue; }

				// Check if part of path (ft. tricky bit flip)
				if (p1.first[!pos] == p2.first[!pos]) {
					path_map[p1.first] = 0;
					break;
				}
			}
		}
	}
}


// DBSCAN Clustering Function, inspired by https://github.com/Eleobert/dbscan/blob/master/dbscan.cpp
std::vector<int> dbscan(ClusterNode *curr_node, const int &points, const int &min_counts,
                        std::vector<std::vector<int>> &assignment, const bool &five) {


	//  BNJ: 7/4/2025 - This function EATS time... I should cache the bounds for each point, make it obvious where to search.
	//						That being said, would only work for the 5' end cause it's ordered.
	//						Could also address the resizing vector issue... 

	int index;
	int dist;
	int clust_num = 0;
	std::vector<int> *adj_vec;
	std::vector<int> neighbors;
	std::vector<int> sub_neighbors;
	std::vector<int> assign_vec(points, -1);
	std::vector<bool> visted(points, false);
	int epsilon = ImpaqtArguments::Args.epsilon;

	// Sepicfy 5' or 3' clusters
	adj_vec = curr_node -> get_five_ref();
	if (!five) { adj_vec = curr_node -> get_three_ref(); }

	// iterate through every point
	for (int i = 0; i < points; i++) {

		neighbors.clear();

		if (visted[i] == false) {

			// Get distance to all other points
			for (int j = 0; j < points; j++) {
				dist = std::abs((*adj_vec)[j] - (*adj_vec)[i]); // Distance between points
				if ((i != j) && (dist <= epsilon)) {
					neighbors.push_back(j);
				}
			}

			// If core point
			if (neighbors.size() >= min_counts) {

				visted[i] = true;
				assign_vec.at(i) = clust_num;
				std::vector<int> cluster_indexes = {i};

				// Continously add points to "stack" and determine if they are core points
				while (neighbors.empty() == false) {

					index = neighbors.back();
					sub_neighbors.clear();
					neighbors.pop_back();

					if (visted[index] == false) {

						visted[index] = true;

						// Check if member of cluster
						for (int k = 0; k < points; k++) {
							if (std::abs((*adj_vec)[index] - (*adj_vec)[k]) <= epsilon) {
								sub_neighbors.push_back(k);
								assign_vec.at(k) = clust_num;
							}
						}

						// If also a core point, copy subneighbors into neighbors to also be checked
						if (sub_neighbors.size() >= min_counts) {
							std::copy(sub_neighbors.begin(), sub_neighbors.end(), std::back_inserter(neighbors));
						}
						cluster_indexes.push_back(index);
					}
				}

				assignment.emplace_back(std::move(cluster_indexes));
				clust_num += 1;
			}
		}
	}
	return assign_vec;
}

// Initiate Transcript Identifying Procedure
void identify_transcripts_dbscan(ClusterList *cluster,  const int &strand) {

	float density;
	bool prime_5 = true;
	int expr, points, min_counts;
	int count_threshold = std::max(ImpaqtArguments::Args.min_count, 10);

	ClusterNode *curr_node = cluster -> get_head(strand);

	while (curr_node != NULL) {

		expr = curr_node -> get_read_count();
		points = curr_node -> get_vec_count();

		// If threshold for transcript detection is reached
		if (expr >= count_threshold) {

			curr_node -> point_sort_vectors();

			std::map<std::string, int> paths;
			std::vector<int> counts;
			std::vector<std::vector<int>> transcripts;
			std::vector<int> assign_vec_5, assign_vec_3;
			std::vector<std::vector<int>> assignments_5,  assignments_3;

			// BNJ - 6/16/2025: Pleaseeeeee fix this casting
			density = (float)expr / (float)(curr_node -> get_stop() - curr_node -> get_start());
			min_counts = std::max((int)((float)expr * (((float)ImpaqtArguments::Args.count_percentage / 100.0))), 10);

			// If not in quantification mode
			if (density < 1.5 || ImpaqtArguments::Args.annotation_file == "") {
				assign_vec_5 = dbscan(curr_node, points, min_counts, assignments_5, prime_5);
				assign_vec_3 = dbscan(curr_node, points, min_counts, assignments_3, !prime_5);

			} else {
				// If read mitochrondrial genome detected, don't bother lol
				curr_node = curr_node -> get_next();
				continue;	
			}
			

			// If clusters were not found
			if (assignments_5.empty() && assignments_3.empty()) {
				;
			
			} else {

				// If clusters were  found
				get_linked_clusters(curr_node, paths, assign_vec_5, assign_vec_3);

				get_coordinates(curr_node, paths,
				                assignments_5, assignments_3,
				                &transcripts, &counts);

                overlap_clusters(curr_node, transcripts, counts);

				report_transcripts(curr_node, transcripts, counts);
				
			}
		}
		curr_node = curr_node -> get_next();
	}
}


void merge_transcripts(ClusterNode *c_node, ClusterNode *n_node) {

	std::vector<std::vector<int>> transcripts = *(c_node -> get_transcripts());
	std::vector<int> counts(c_node -> get_transcript_num(), 0);

	// Populate New Transcript and Count Vecs
	for (int i = 0; i < c_node -> get_transcript_num(); i++) { counts[i] = (int)(c_node -> get_transcript_expr(i)); }
	for (int i = 0; i < n_node -> get_transcript_num(); i++) {
		transcripts.push_back(n_node -> get_transcripts() -> at(i));
		counts.push_back((int)(n_node -> get_transcript_expr(i)));
	}
	c_node -> clear_transcripts();

	// Merge Final Transcripts
	overlap_clusters(c_node, transcripts, counts);
	report_transcripts(c_node, transcripts, counts);
}


// Combine clusters with nonzero neighbors
void collapse_transcripts(ClusterList *cluster, int t_strand) {

	ClusterNode *c_node = cluster -> get_head(t_strand);
	ClusterNode *n_node = NULL;

	while (c_node != NULL) {

		// Only Collapse transcripts
		if (c_node -> get_transcript_num() != 0) {

			n_node = c_node -> get_next();

			while (n_node != NULL) {

				if (n_node -> get_transcript_num() != 0) {

					if (c_node -> get_transcript_stop() >= n_node -> get_transcript_start()) {

						merge_transcripts(c_node, n_node);
						c_node -> update_read_counts(n_node -> get_read_count());
						c_node -> update_vec_counts(n_node -> get_vec_count());
						c_node -> update_stop(n_node -> get_stop());
						n_node -> set_skip();

					} else { break; }
				}
				
				n_node = n_node -> get_next();
			}
			c_node -> quantify_transcripts();

		} else { n_node = c_node -> get_next(); }

		c_node = n_node;
	}
}

void identify_transcripts(ClusterList *cluster, int t_strand) {
	identify_transcripts_dbscan(cluster, t_strand);
	collapse_transcripts(cluster, t_strand);
}