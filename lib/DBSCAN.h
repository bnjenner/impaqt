//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DBSCAN and Related Functions

// Check if transcripts overlap / are contained in another transcript
bool check_subset(const std::vector<int>& a, const std::vector<int>& b) {

	bool match = false;
	bool started = false;

	int i = 0;
	int j = 0;
	int n = a.size() / 2;
	int m = b.size() / 2;

	// Iterate through exons, see if match is complete 
	while (j < m) {

		// If end of A is reached or first match was not a hit
		if (i >= n || (!match && ((i > 0) && (j > 0)))) {
			break;

		} else if (check_bounds(a[(i*2)], a[(i*2)+1], b[(j*2)], b[(j*2)+1])) {
			match = true; started = true;
			i += 1;

			// If last exon and it matches the first exon of B is after but close to A
		} else if (i == n - 1 && (a.back() < b[0]) && 
				   std::abs(a.back() - b[0]) <= ImpaqtArguments::Args.epsilon) {
			match = true; started = true;
			i += 1;

		} else {

			if (started) {
				match = false; break;
			} else {
				i += 1; j -= 1;
			}
		}
		j += 1;
	}

	return match;
}


// Get Core Points of Transcript
int get_quant(const std::vector <int> result, const std::vector<std::vector<int>> &init_copy, const std::vector<int> counts) {
	int core_points = 0;
	for (int i = 0; i < init_copy.size(); i++) {
		if (check_subset(result, init_copy[i])) {
			core_points += counts[i];
		}
	}
	return core_points;
}


// Merge Overlapping Transcripts
std::vector<int> merge_transcripts(const std::vector<int>& a, const std::vector<int>& b) {

	std::vector<int> merged;

	int i = 0;
	int j = 0;
	int n = a.size() / 2;
	int m = b.size() / 2;

	bool appended = false;

	while (j < m) {

		// If section of A precedes B
		if (i != n && j == 0 && (a[(i*2)+1] < b[(j*2)])) {
			merged.push_back(a[(i*2)]);
			merged.push_back(a[(i*2)+1]);
			i += 1;

			// If remaining sections in B, but A is spent
		} else if (i == n) {

			if (!appended && std::abs(merged.back() - b[(j*2)]) <= ImpaqtArguments::Args.epsilon) {
				merged.back() = std::max(merged.back(), b[(j*2)+1]);
				appended = true;
			} else {
				merged.push_back(b[(j*2)]);
				merged.push_back(b[(j*2)+1]);
			}
			j += 1;
		
			// Merge Overlapping sections
		} else {
			merged.push_back(std::min(a[(i*2)], b[(j*2)]));
			merged.push_back(std::max(a[(i*2)+1], b[(j*2)+1]));
			i += 1;
			j += 1;
		}
	}

	// If remaining sections of A
	while (i < n) {
		merged.push_back(a[(i*2)]);
		merged.push_back(a[(i*2)+1]);
		i += 1;
	}
	
	return merged;
}

// Reduce Transcript Number by Overlapping. Report Unique Transcripts
void get_final_transcripts(ClusterNode *curr_node, std::vector<std::vector<int>> &transcripts, std::vector<int> &counts) {

	int n;
	bool all_unique, overlap;
	std::vector<int> absorbed;

	std::vector<int> tmp;
	std::vector<std::vector<int>> result;
	std::vector<std::vector<int>> init_copy = transcripts;

	// No need to overlap
	if (transcripts.size() == 1) { 
		curr_node -> add_transcript(transcripts[0], counts[0]);
		return;
	}

	// Reverse and Negative if reverse strand
	if (curr_node -> get_strand() == 1) { reverse_transcripts(transcripts); }

	// Merge transcripts until all are unique
	while (true) {

		// Reset Checks
		n = transcripts.size();
		result.clear(); 
		absorbed.clear();
		all_unique = true;

		// Iterate through transcripts
		for (int i = 0; i < n; i++) {
			overlap = false;

			for (int j = i + 1; j < n; j++) {

				// Check if transcripts overlap
				if (check_subset(transcripts[i], transcripts[j])) {

					overlap = true;
					absorbed.push_back(j);
					tmp = merge_transcripts(transcripts[i], transcripts[j]);
					
					// If new transcript not already in results, add
					auto it = std::find(result.begin(), result.end(), tmp);
					if (it == result.end()) { result.push_back(tmp); }
				}
			}

			// If unique transcripts and not absorbed, add. If not, add another iteration
			if (!overlap && std::find(absorbed.begin(), absorbed.end(), i) == absorbed.end()) {
				result.push_back(transcripts[i]);
			} else {
				all_unique = false;
			}
		}

		if (!all_unique) {
			transcripts = result;
		} else {
			break;
		}
	}


	// Reverse and Negative Results if Necessary
	if (curr_node -> get_strand() == 1) { reverse_transcripts(result); }
	
	// Report Final Transcripts
	int core_points = 0;
	for (int i = 0; i < result.size(); i++) {
		core_points = get_quant(result[i], init_copy, counts);
		curr_node -> add_transcript(result[i], core_points);
	}
}


// Report Unique Transcripts (no overlapping, used for Mitochrondria)
void report_transcripts(ClusterNode *curr_node, std::vector<std::vector<int>> &result, std::vector<int> &counts) {
	
	// Reverse and Negative Results if Necessary
	if (curr_node -> get_strand() == 1) { reverse_transcripts(result); }

	// Report Final Transcripts
	for (int i = 0; i < result.size(); i++) {
		curr_node -> add_transcript(result[i], counts[i]);
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

		if (p.second == 0) { continue; } // Skip if no reads

		std::vector<int> tmp_vec;

		// add 5' region
		if (p.first[0] != '-') {

			// Get Bounds of cluster
			index = std::stoi(p.first.substr(0, 1));
			min_pos = get_pos_min(index, core_5, curr_node -> get_five_ref());
			max_pos = get_pos_max(index, core_5, curr_node -> get_five_ref());

			// swap variables if necessary
			if (min_pos > max_pos) { variable_swap(min_pos, max_pos); } 

			tmp_vec.push_back(min_pos);
			tmp_vec.push_back(max_pos);
		}

		// add 3' region
		if (p.first[1] != '-') {

			// Get Bounds of cluster
			index = std::stoi(p.first.substr(1, 1));
			min_pos = get_pos_min(index, core_3, curr_node -> get_three_ref());
			max_pos = get_pos_max(index, core_3, curr_node -> get_three_ref());

			// swap variables if necessary
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

		transcripts -> push_back(tmp_vec);
		counts -> push_back(p.second);
	}
}


// Find all linked DBSCAN clusters
void get_linked_clusters(ClusterNode *curr_node, std::map<std::string, int> &path_map,
                         const std::vector<int> &assign_5, const std::vector<int> &assign_3) {

	std::string path;

	// Use hashmap to store clusters and supporting counts
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

		// Add to map and increment
		if (path_map.find(path) != path_map.end()) {
			path_map[path] += 1; // Increment count if path already exists
		} else {
			path_map[path] = 1; // Add new path with count of 1
		}
	}


	// Absorb orphan paths
	int pos;
	for (const auto& p1 : path_map) {

		if (path_map[p1.first] == 0) { continue; }

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


// DBSCAN Clustering Function
// 		inspired by https://github.com/Eleobert/dbscan/blob/master/dbscan.cpp
std::vector<int> dbscan(ClusterNode *curr_node, const int &points, const int &min_counts,
                        std::vector<std::vector<int>> &assignment, const bool &five, const bool &mito) {

	int index;
	int dist;
	int clust_num = 0;
	std::vector<int> *adj_vec;
	std::vector<int> neighbors;
	std::vector<int> sub_neighbors;
	std::vector<int> assign_vec(points, -1);
	std::vector<bool> visted(points, false);


	int epsilon = ImpaqtArguments::Args.epsilon;
	if (mito) { epsilon = 50; } // If mito, use smaller epsilon (magic number again... look they're fundamentally different problems)

	// Sepicfy 5' or 3' clusters
	if (five) {
		adj_vec = curr_node -> get_five_ref();
	} else {
		adj_vec = curr_node -> get_three_ref();
	}


	// iterate through every point
	for (int i = 0; i < points; i++) {

		neighbors.clear();

		if (visted[i] == false) {

			// Get distance to all other points
			for (int j = 0; j < points; j++) {
				dist = std::abs((*adj_vec)[j] - (*adj_vec)[i]); // Distance between points
				if ((i != j) && (dist <= epsilon)) { neighbors.push_back(j); }
			}

			// If core point
			if (neighbors.size() >= min_counts) {

				visted[i] = true;
				assign_vec.at(i) = clust_num;
				std::vector<int> cluster_indexes;

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
void find_transcripts_DBSCAN(ClusterList &cluster,  const int &strand) {

	float density;
	int expr, points, min_counts;
	int count_threshold = std::max(ImpaqtArguments::Args.min_count, 10);

	ClusterNode *curr_node = cluster.get_head(strand);

	while (curr_node != NULL) {

		expr = curr_node -> get_read_count();
		points = curr_node -> get_vec_count();

		// If threshold for transcript detection is reached
		if (expr >= count_threshold) {

			std::map<std::string, int> paths;
			std::vector<int> counts;
			std::vector<std::vector<int>> transcripts;
			std::vector<int> assign_vec_5, assign_vec_3;
			std::vector<std::vector<int>> assignments_5,  assignments_3;

			// Read Density and Min Counts for DBSCAN
			density = (float)expr / (float)(curr_node -> get_stop()) - curr_node -> get_start();
			min_counts = std::max((int)((float)expr * (((float)ImpaqtArguments::Args.count_percentage / 100.0))), 10);

			// Sort Vectors
			curr_node -> sort_vectors();

			// Run DBSCAN
			if (density < 1.5) {
				assign_vec_5 = dbscan(curr_node, points, min_counts, assignments_5, true, false);
				assign_vec_3 = dbscan(curr_node, points, min_counts, assignments_3, false, false);

			} else {

				// If read mitochrondrial genome detected, only cluster 5' end 
				min_counts = (int)((float)expr * 0.01); // 1% of total reads for mito (magic number, I am sorry)
				assign_vec_5 = dbscan(curr_node, points, min_counts, assignments_5, true, true);
				assign_vec_3 = std::vector<int>(points, -1);
				assignments_3.push_back(assign_vec_3);		
			}
			

			// If clusters were not found
			if (assignments_5.empty() && assignments_3.empty()) {
				;
			
			} else {

				// Find all linked DBSCAN clusters
				get_linked_clusters(curr_node, paths, assign_vec_5, assign_vec_3);

				// Get Transcript Coords and Core Points
				get_coordinates(curr_node, paths,
				                assignments_5, assignments_3,
				                &transcripts, &counts);

				
				// Report Final Transcripts
				if (density < 1.5) {
					// Merge overlapping transcripts 
					get_final_transcripts(curr_node, transcripts, counts);
				} else {
					// Handle Mitochrondria, do not merge
					report_transcripts(curr_node, transcripts, counts);
				}

				// Determine abundance of each transcript
				curr_node -> quantify_transcripts();
			}

		}

		curr_node = curr_node -> get_next();
	}
}