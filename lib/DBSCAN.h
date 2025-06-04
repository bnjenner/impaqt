//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// DBSCAN and Related Functions

// Reverse and Negate Reverse Strand Vectors
std::vector<int> reverse_and_negate(const std::vector<int> &vec) {
	int x = 0;
	int n = vec.size();
	std::vector<int> tmp_vec(n);
	for (int j = n - 1; j > -1; j--) {
		tmp_vec[x] = vec[j] * -1;
		++x;
	}
	return tmp_vec;
}

// Used to Sort those Pesky Reverse Strands
bool compare_first_element(const std::vector<int>& a, const std::vector<int>& b) { return a[0] < b[0]; }

// Sort Using Vector Lengths
bool compare_length(const std::vector<int>& a, const std::vector<int>& b) { return a.size() < b.size(); }

// Checks to see if regions of transcript overlap
bool check_bounds(const int &a_start, const int &a_stop, const int &b_start, const int &b_stop) {
	// Precedes
	if (a_stop >= b_start && a_stop <= b_stop) {
		return true;
		// Spans
	} else if (a_start <= b_start && a_stop >= b_stop) {
		return true;
		// Follows
	} else if (a_start >= b_start && a_start <= b_stop) {
		return true;
		// Precedes (but close enough)
	} else if (std::abs(a_stop - b_start) <= ImpaqtArguments::Args.epsilon) {
		return true;
	}
	return false;
}

// Check if Transcript is a subset of another transcript
bool check_if_subset(const std::vector<int>& a, const std::vector<int>& b) {

	int a_len = a.size();
	int b_len = b.size();
	if (a_len > b_len) { return false; }

	bool started = false;
	bool match = false;

	for (int i = 0; i < (a_len / 2); i++) {
		for (int j = 0; j < (b_len / 2); j++) {

			if (check_bounds(a[(2 * i)], a[(2 * i) + 1], b[(2 * j)], b[(2 * j) + 1])) {
				if (!started) { started = true; } // If path tracing is started
				match = true;
				++i;
			} else {
				match = false;
				if (started) { break; }
			}

			if (i >= (a_len) / 2) {
				break;
			}
		}

		// Only return true if a full match is achieved
		if (started && match) {
			return true;
		} else {
			break;
		}
	}
	return false;
}

// Merge Overlapping Transcripts
std::vector<int> reduce_transcripts(const std::vector<int>& a, const std::vector<int>& b) {

	int a_len = a.size();
	int b_len = b.size();
	std::vector<int> result(a);

	for (int i = 0; i < (b_len / 2); i++) {
		if (b[(2 * i)] > a[a_len - 1]) {
			result.push_back(b[(2 * i)]);
			result.push_back(b[(2 * i) + 1]);
		} else {
			if (b[(2 * i)] < a[a_len - 2]) { result[a_len - 2] = b[(2 * i)]; }
			if (b[(2 * i) + 1] > a[a_len - 1]) { result[a_len - 1] = b[(2 * i) + 1]; }
		}
	}
	return result;
}


// Merge Overlapping Transcripts
void get_final_transcripts(ClusterNode *curr_node, std::vector<std::vector<int>> &transcripts, std::vector<float> &counts) {

	// BNJ: 5/31/2025 - I can probably find a better way to do this lol

	int n = transcripts.size();
	std::vector<std::vector<int>> result;
	std::vector<float> res_counts;

	// Reverse and Negative if reverse strand (I'm actually pretty proud of this solution)
	if (curr_node -> get_strand() == 1) {
		for (int i = 0; i < transcripts.size(); i++) {
			transcripts[i] = reverse_and_negate(transcripts[i]);
		}
		std::sort(transcripts.begin(), transcripts.end(), compare_first_element);
	}

	// Overlap transcripts
	int t_len;
	bool overlap;
	for (int i = 0; i < n; i++) {
		overlap = false;
		t_len = transcripts[i].size();

		for (int j = i + 1; j < n; j++) {

			if (check_bounds(transcripts[i][t_len - 2], transcripts[i][t_len - 1], transcripts[j][0], transcripts[j][1])) {
				result.push_back(reduce_transcripts(transcripts[i], transcripts[j])); // Merge transcript and add to results
				res_counts.push_back(counts[i] + counts[j]);
				overlap = true;
			}
		}
		if (!overlap) {
			result.push_back(transcripts[i]);
			res_counts.push_back(counts[i]);
		}
	}

	bool unique;
	bool all_unique = false;

	// Further Clean Up, some transcripts are subsets of others, iterate until all gone
	while (!all_unique) {

		// Sort by size and create tmp vec
		all_unique = true;
		n = result.size();
		std::sort(result.begin(), result.end(), compare_length);
		std::vector<std::vector<int>> tmp_result;
		std::vector<float> tmp_counts;

		// Check for transcripts that subset another
		for (int i = 0; i < n; i++) {
			unique = true;
			for (int j = i + 1; j < n; j++) {
				if (check_if_subset(result[i], result[j])) {
					all_unique = false;
					unique = false;
					break;
				}
			}

			if (unique) {
				tmp_result.push_back(result[i]);
				tmp_counts.push_back(res_counts[i]);
			}
		}
		// If merges made
		if (!all_unique) {
			result = tmp_result;
			res_counts = tmp_counts;
		}
	}

	// Report Final Transcripts
	for (int i = 0; i < result.size(); i++) {
		// Reinstate original order if necessary
		if (curr_node -> get_strand() == 1) {
			result[i] = reverse_and_negate(result[i]);
		}
		curr_node -> add_transcript(result[i], res_counts[i]);
	}
}

// Get Transcript Coordinates
void get_coordinates(ClusterNode *curr_node, const std::vector<std::string> &paths,
                     std::vector<std::vector<int>> &core_5, std::vector<std::vector<int>> &core_3,
                     std::vector<std::vector<int>> *transcripts, std::vector<float> *counts) {

	float core_points;
	int index, min_pos, max_pos;  // necessary because 3' vec may not be in order due to splicing

	// Get transcript coordinates
	//	BNJ: 5/2/2025 - Horrible extending vector, but they're short so it shouldn't matter THAT much
	//	BNJ: 5/31/2025 - Also worth mentioning, the tmp_Vec should never be more than 4 in length
	//	BNJ: 6/4/2025 - Man, I should really use bits instead of a string here. Or maybe not?

	for (const auto &p : paths) {

		core_points = 0.0f;
		std::vector<int> tmp_vec;
		std::vector<int>::iterator min_result;
		std::vector<int>::iterator max_result;

		// add 5' region
		if (p[0] != '-') {

			// Get Bounds of cluster
			index = std::stoi(p.substr(0, 1));
			min_result = std::min_element(core_5.at(index).begin(), core_5.at(index).end());
			max_result = std::max_element(core_5.at(index).begin(), core_5.at(index).end());
			min_pos = curr_node -> get_five_vec().at(*min_result);
			max_pos = curr_node -> get_five_vec().at(*max_result);

			// swap variables if necessary
			if (min_pos > max_pos) {
				max_pos = max_pos + min_pos;
				min_pos = max_pos - min_pos;
				max_pos = max_pos - min_pos;
			}

			tmp_vec.push_back(min_pos);
			tmp_vec.push_back(max_pos);
			core_points += (float)core_5.at(index).size();

		}

		// add 3' region
		if (p[1] != '-') {

			// Get Bounds of cluster
			index = std::stoi(p.substr(1, 1));
			min_result = std::min_element(core_3.at(index).begin(), core_3.at(index).end());
			max_result = std::max_element(core_3.at(index).begin(), core_3.at(index).end());
			min_pos = curr_node -> get_three_vec().at(*min_result);
			max_pos = curr_node -> get_three_vec().at(*max_result);

			// swap variables if necessary
			if (min_pos > max_pos) {
				max_pos = max_pos + min_pos;
				min_pos = max_pos - min_pos;
				max_pos = max_pos - min_pos;
			}

			tmp_vec.push_back(min_pos);
			tmp_vec.push_back(max_pos);
			core_points += (float)core_3.at(index).size();

		}

		// Reorder if necessary (really just a reverse strand thing, will figure this out)
		if (curr_node -> get_strand() == 1 && tmp_vec.size() > 2
		        && tmp_vec[0] >= tmp_vec[2]) {

			// If completely encompassing
			if (tmp_vec[1] <= tmp_vec[3]) {
				tmp_vec = {tmp_vec[2], tmp_vec[3]};

				// if separate clusters
			} else if (tmp_vec[3] < tmp_vec[0]) {
				tmp_vec = {tmp_vec[2], tmp_vec[3],
				           tmp_vec[0], tmp_vec[1]
				          };

				// if overlapping
			} else {
				tmp_vec = {tmp_vec[2], tmp_vec[1]};
			}
		}

		// if close enough
		if (tmp_vec.size() > 2 && ImpaqtArguments::Args.epsilon >= (tmp_vec[2] - tmp_vec[1])) {
			tmp_vec = {tmp_vec[0], tmp_vec[3]};
		}

		transcripts -> push_back(tmp_vec);
		counts -> push_back(core_points);
	}
}


// Find all linked DBSCAN clusters
void get_linked_clusters(ClusterNode *curr_node, std::vector<std::string> &path_vec,
                         const std::vector<int> &assign_5, const std::vector<int> &assign_3) {

	std::string path;
	std::vector<std::string> tmp_vec;

	for (int i = 0; i < curr_node -> get_read_count(); i++) {
		path = "";

		// assigned in 5' DBSCAN
		if (assign_5.at(i) != -1) {

			path = std::to_string(assign_5.at(i)) + '-';
			if (assign_3.at(i) != -1) { path.at(1) = std::to_string(assign_3.at(i))[0]; }

			// unassigned in 5' DBSCAN
		} else if (assign_5.at(i) == -1 && assign_3.at(i) != -1) {
			path = "-" + std::to_string(assign_3.at(i));
		}

		// add to paths vector if path is not empty
		auto it = std::find(tmp_vec.begin(), tmp_vec.end(), path);
		if (path != "" && it == tmp_vec.end()) { tmp_vec.push_back(path); }
	}


	// Remove orphan paths that are part of other paths
	int pos;
	bool add;
	for (int i = 0; i < tmp_vec.size(); i++) {

		add = true;
		pos = tmp_vec[i].find('-');

		// if path is orphaned
		if (pos != std::string::npos) {
			for (int j = 0; j < tmp_vec.size(); j++) {

				if (i == j) { continue; }

				// Check if part of path (ft. tricky bit flip)
				if (tmp_vec[i][!pos] == tmp_vec[j][!pos]) {
					add = false;
					break;
				}
			}
		}

		if (add) { path_vec.push_back(tmp_vec[i]); }
	}
}


// DBSCAN Clustering Function
// 		inspired by https://github.com/Eleobert/dbscan/blob/master/dbscan.cpp
std::vector<int> dbscan(ClusterNode *curr_node, const int &points, const int &min_counts,
                        std::vector<std::vector<int>> &assignment, const bool &five) {

	int index;
	int dist;
	int clust_num = 0;
	std::vector<int> *adj_vec;
	std::vector<int> neighbors;
	std::vector<int> sub_neighbors;
	std::vector<int> assign_vec(points, -1);
	std::vector<bool> visted(points, false);


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
				if ((i != j) && (dist <= ImpaqtArguments::Args.epsilon)) { neighbors.push_back(j); }
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
							if (std::abs((*adj_vec)[index] - (*adj_vec)[k]) <= ImpaqtArguments::Args.epsilon) {
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

// Initiate DBSCAN Transcript Identifying Procedure
void find_transcripts_DBSCAN(ClusterList &cluster,  const int &strand) {

	int points;
	int min_counts;
	int count_threshold = std::max(ImpaqtArguments::Args.min_count, 10);

	ClusterNode *curr_node = cluster.get_head(strand);

	while (curr_node != NULL) {

		points = curr_node -> get_read_count();

		// If threshold for transcript detection is reached
		if (points >= count_threshold) {

			std::vector<std::string> paths;
			std::vector<float> counts;
			std::vector<std::vector<int>> transcripts;
			std::vector<int> assign_vec_5, assign_vec_3;
			std::vector<std::vector<int>> assignments_5,  assignments_3;

			// Min Counts for DBSCAN
			min_counts = std::max((int)((float)points * ((float)(ImpaqtArguments::Args.count_percentage / 100))), 20);

			// Run DBSCAN
			assign_vec_5 = dbscan(curr_node, points, min_counts, assignments_5, true);
			assign_vec_3 = dbscan(curr_node, points, min_counts, assignments_3, false);


			// If clusters were not found
			if (assignments_5.empty() && assignments_3.empty()) {
				curr_node = curr_node -> get_next();
				continue;
			}

			// Find all linked DBSCAN clusters
			get_linked_clusters(curr_node, paths, assign_vec_5, assign_vec_3);

			// Get Transcript Coords and Core Points
			get_coordinates(curr_node, paths,
			                assignments_5, assignments_3,
			                &transcripts, &counts);

			// Merge overlapping transcripts and add to cluster node
			get_final_transcripts(curr_node, transcripts, counts);
		}

		curr_node = curr_node -> get_next();
	}
}