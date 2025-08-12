#include <iostream>
#include <iomanip>

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
void overlap_clusters(ClusterNode *node, std::vector<std::vector<int>> &transcripts, std::vector<int> &counts) {

	std::vector<std::vector<int>> init_copy = transcripts;

	// Return if no need to overlap
	if (transcripts.size() == 1) { return; }
	if (node -> get_strand() == 1) { reverse_transcripts(transcripts); } 

	// Overlap Transcripts until unique
	bool unique = overlap_aux(transcripts);
	while (!unique) {
		unique = overlap_aux(transcripts);
	}

	// Reverse and Negative Results if Necessary
	if (node -> get_strand() == 1) { reverse_transcripts(transcripts); } 
	std::sort(transcripts.begin(), transcripts.end(), compare_first_element);
	
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
void get_coordinates(ClusterNode *node, const std::map<std::string, int> &paths,
                     const std::map<int, std::vector<int>> &regions_5, const std::map<int, std::vector<int>> &regions_3,
                     std::vector<std::vector<int>> *transcripts, std::vector<int> *counts) {

	// Get transcript coordinates
	//	BNJ: 5/31/2025 - Worth mentioning, the tmp_vec should never be more than 4 in length
	//	BNJ: 6/4/2025 - Man, I should really find a better struct than strings for the paths.
	//  BNJ: 8/10/2025 - ... Say that again? 

	int index, n; 

	for (const auto &p : paths) {

		if (p.second == 0) { continue; }

		std::vector<int> tmp_vec;

		// add 5' region
		if (p.first[0] != '-') {
			index = std::stoi(p.first.substr(0, 1));
			tmp_vec.push_back(regions_5.at(index)[0]);
			tmp_vec.push_back(regions_5.at(index)[1]);
			if (tmp_vec[0] > tmp_vec[1]) { variable_swap(tmp_vec[0], tmp_vec[1]); } 
		}

		// add 3' region
		if (p.first[1] != '-') {
			n = tmp_vec.size() + 2;
			index = std::stoi(p.first.substr(1, 1));
			tmp_vec.push_back(regions_3.at(index)[0]);
			tmp_vec.push_back(regions_3.at(index)[1]);
			if (tmp_vec[n-2] > tmp_vec[n-1]) { variable_swap(tmp_vec[n-2], tmp_vec[n-1]); } 
		}

		if (tmp_vec.size() > 2) {
			// If out of order
			if (tmp_vec[2] < tmp_vec[0]) {
				tmp_vec = {tmp_vec[2], tmp_vec[3], tmp_vec[0], tmp_vec[1]};
			}
			// If two regions and they are close or out of order, merge
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
void get_linked_clusters(ClusterNode *node, std::map<std::string, int> &path_map,
                         const std::vector<int> &assign_5, const std::vector<int> &assign_3) {

	std::string path;

	for (int i = 0; i < node -> get_vec_count(); i++) {

		path = "";

		// assigned in 5' DBSCAN
		if (assign_5.at(i) != -1) {
			path = std::to_string(assign_5.at(i)) + '-';
			if (assign_3.at(i) != -1) {
				path.at(1) = std::to_string(assign_3.at(i))[0];
			}

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

		if (p1.second < 10) { path_map[p1.first] = 0; continue; }

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


// Get Nearest Neighbors 
std::vector<int> get_nearest_neighbors(const int &i, const int &points, std::vector<bool> &queued,
                           			   const std::vector<int> &indices, const std::vector<int> *adj_vec) {

	int bound;
	int p = indices[i];
	std::vector<int> neighbors;

	// Forward Search
	bound = (adj_vec -> at(p)) + ImpaqtArguments::Args.epsilon;
	for (int j = i + 1; j < points; j++) {
		if (adj_vec -> at(indices[j]) > bound) { break; }
		if (!queued.empty()) { queued[j] = true; }
		neighbors.push_back(j);
	}

	// Backward Search
	bound = (adj_vec -> at(p)) - ImpaqtArguments::Args.epsilon;
	for (int j = i - 1; j >= 0; j--) {
		if (adj_vec -> at(indices[j]) < bound) { break; }
		if (!queued.empty()) { queued[j] = true; }
		neighbors.push_back(j);
	}

	return neighbors;
} 


// DBSCAN Clustering Function, inspired by https://github.com/Eleobert/dbscan/blob/master/dbscan.cpp
std::vector<int> dbscan(ClusterNode *node, const int &points, const int &min_counts,
                        std::map<int, std::vector<int>> &regions, const bool &five) {


	// DBSCAN Variables
 	bool skip;
	int clust_num = 0;
	int p1, p2, index;
	int min_point, max_point;
	std::vector<int> *adj_vec;
	std::vector<bool> empty_vec;
	std::vector<int> neighbors, sub_neighbors;

	// Indexing Variables
	std::vector<int> indices(points);
	std::vector<int> assign_vec(points, -1);
	std::vector<bool> visted(points, false);


	adj_vec = node -> get_five_ref();
	std::iota(indices.begin(), indices.end(), 0);

	// If dealing with the unsorted vector, get sorted indices
	if (!five) {
		adj_vec = node -> get_three_ref();
		std::sort(indices.begin(), indices.end(),
			       [&](int i, int j) -> bool {
			            return (*adj_vec)[i] < (*adj_vec)[j];
			        }
		         );
	}


	for (int i = 0; i < points; i++) {

		if (visted[i] == true) { continue; }

		p1 = indices[i];
		std::vector<bool> queued(points, false);
		
		neighbors = get_nearest_neighbors(i, points, queued, indices, adj_vec);

		// If core point
		if (neighbors.size() >= min_counts) {

			visted[i] = true;
			assign_vec.at(p1) = clust_num;
			std::vector<int> cluster_points = {adj_vec -> at(p1)};

			int x = 0;
			while (x < neighbors.size()) {

				index = neighbors[x];
				p2 = indices[index];

				if (visted[index] == false) {

					// Skip Duplicate Points
					skip = false;
					for (const auto &t_point : cluster_points) {
						if (adj_vec -> at(p2) == t_point) {
							assign_vec.at(p2) = clust_num;
							visted[index] = true;
							skip = true; 
							break;
						}
					}

					if (!skip) {

						visted[index] = true;
						assign_vec.at(p2) = clust_num;

						sub_neighbors = get_nearest_neighbors(index, points, empty_vec, indices, adj_vec);

						// If also a core point, copy subneighbors into neighbors to also be checked
						if (sub_neighbors.size() >= min_counts) {
							for (const auto &n : sub_neighbors) {
								if (!queued[n]) {
									neighbors.push_back(n);
									queued[n] = true;
								}
							}
							cluster_points.push_back(adj_vec -> at(p2));
						}
					}					
				}
				++x;			
			}

			min_point = *std::min_element(cluster_points.begin(), cluster_points.end());
			max_point = *std::max_element(cluster_points.begin(), cluster_points.end());
			regions[clust_num] = std::vector<int>{min_point, max_point};
			clust_num += 1;
			i = x - 1; // Skip to next unvisited point
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

	std::map<std::string, int> paths;
	std::vector<int> counts;
	std::vector<std::vector<int>> transcripts;
	std::vector<int> assign_vec_5, assign_vec_3;
	std::map<int, std::vector<int>> regions_5, regions_3;

	ClusterNode *node = cluster -> get_head(strand);

	while (node != NULL) {

		expr = node -> get_read_count();
		points = node -> get_vec_count();

		// If threshold for transcript detection is reached
		if (expr >= count_threshold) {

			// Reset Data
			paths.clear();
			regions_5.clear();
			regions_3.clear();
			transcripts.clear();
			counts.clear();

			// Sort vectors
			node -> point_sort_vectors();

			// BNJ - 6/16/2025: Pleaseeeeee fix this casting
			density = (float)expr / (float)(node -> get_stop() - node -> get_start());
			min_counts = std::max((int)((float)expr * (((float)ImpaqtArguments::Args.count_percentage / 100.0))), 10);

			// If not in quantification mode
			if (density < ImpaqtArguments::Args.density_threshold || ImpaqtArguments::Args.annotation_file == "") {
				assign_vec_5 = dbscan(node, points, min_counts, regions_5, prime_5);
				assign_vec_3 = dbscan(node, points, min_counts, regions_3, !prime_5);

			} else {		
				// If read mitochrondrial genome detected, don't bother lol
				std::cerr << "//    NOTICE: Density threshold met (" 
						  << std::fixed << std::setprecision(2) << density << "). Skipping "
						  << node -> get_contig_name() << ":" 
						  << node -> get_start() << "-" 
						  << node -> get_stop() << ".\n";
				node = node -> get_next();
				continue;	
			}
			
			// If clusters were not found
			if (regions_5.empty() && regions_3.empty()) {
				;
			
			} else {

				// If clusters were  found
				get_linked_clusters(node, paths, assign_vec_5, assign_vec_3);

				get_coordinates(node, paths,
				                regions_5, regions_3,
				                &transcripts, &counts);


				// If no transcripts have at least 10 supporting reads. (maybe don't hardcode this?)
				if (!transcripts.empty()) {
					overlap_clusters(node, transcripts, counts);
					report_transcripts(node, transcripts, counts);
					node -> empty_vectors();
				}
			}
		}
		node = node -> get_next();
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