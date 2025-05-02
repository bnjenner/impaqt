#include <algorithm>
#include <string>

void reduce_transcripts(ClusterNode *curr_node, const std::vector<std::string> &paths,
                        std::vector<std::vector<int>> &core_5, std::vector<std::vector<int>> &core_3) {

	std::vector<std::vector<int>> tmp_transcripts;

	// Horrible extending vector, but it's short so it shouldn't matter THAT much
	for (const auto &p : paths) {

		int index;
		int min_pos, max_pos; // necessary because 3' vec may not be in order due to splicing
		std::vector<int> tmp_vec;
		std::vector<int>::iterator min_result;
		std::vector<int>::iterator max_result;

		// add 5' region
		if (p[0] != '-') {
			index = std::stoi(p.substr(0, 1));
			min_result = std::min_element(core_5.at(index).begin(), core_5.at(index).end());
			max_result = std::max_element(core_5.at(index).begin(), core_5.at(index).end());
			tmp_vec.push_back(curr_node -> get_five_vec().at(*min_result));
			tmp_vec.push_back(curr_node -> get_five_vec().at(*max_result));
		}

		// add 3' region
		if (p[1] != '-') {
			index = std::stoi(p.substr(1, 1));
			min_result = std::min_element(core_3.at(index).begin(), core_3.at(index).end());
			max_result = std::max_element(core_3.at(index).begin(), core_3.at(index).end());

			min_pos = curr_node -> get_three_vec().at(*min_result);
			max_pos = curr_node -> get_three_vec().at(*max_result);

			// swap variables
			if (min_pos > max_pos) {
				max_pos = max_pos + min_pos;
				min_pos = max_pos - min_pos;
				max_pos = max_pos - min_pos;
			}

			// if clusters within read length
			if (!tmp_vec.empty() && ImpaqtArguments::Args.epsilon >= (min_pos - tmp_vec.back())) {
				tmp_vec[1] = max_pos;

			} else {
				tmp_vec.push_back(min_pos);
				tmp_vec.push_back(max_pos);
			}
		}

		tmp_transcripts.push_back(tmp_vec);
	}

	// for (int i = 0; i < paths.size(); i++) {
	// 	std::cerr << paths[i] << "\t";
	// 	for (int j = 0; j < tmp_transcripts[i].size(); j++) {
	// 		std::cerr << tmp_transcripts[i][j] << " ";
	// 	}
	// 	std::cerr << "\n";
	// }


	/*
	Ok so there is really only one way clusters can overlap
		1.) The chain, 3' of one cluster overlaps with 5' of following cluster.

	*/

	int pos;
	for (int i = 0; i < tmp_transcripts.size() - 1; i++) {

		for (int j = i + 1; j < tmp_transcripts.size(); j++) {

			pos = tmp_transcripts[i].back();

			if (pos >= tmp_transcripts[j][0] && pos <= tmp_transcripts[j][1]) {

				tmp_transcripts[i].back() = tmp_transcripts[j][1];

				if (tmp_transcripts[j].size() > 2) {
					tmp_transcripts[j].push_back(tmp_transcripts[j][2]);
					tmp_transcripts[j].push_back(tmp_transcripts[j][3]);
				}

				tmp_transcripts[j].clear();
				i = j;
			}

			j += 1;
		}
	}


	// for (int i = 0; i < paths.size(); i++) {
	// 	std::cerr << paths[i] << "\t";

	// 	if (tmp_transcripts[i].empty()) {
	// 		std::cerr << "MERGED\n";
	// 		continue;
	// 	}

	// 	for (int j = 0; j < tmp_transcripts[i].size(); j++) {
	// 		std::cerr << tmp_transcripts[i][j] << " ";
	// 	}
	// 	std::cerr << "\n";
	// }

}

void trace_transcripts(ClusterNode *curr_node, std::vector<std::string> &path_vec,
                       const std::vector<int> &assign_5, const std::vector<int> &assign_3) {

	std::string path;
	std::vector<std::string> tmp_vec;

	for (int i = 0; i < curr_node -> get_read_count(); i++) {
		path = "";

		// assigned in 5' DBSCAN
		if (assign_5.at(i) != -1) {

			// assigned in 3' DBSCAN
			if (assign_3.at(i) != -1) {
				path = std::to_string(assign_5.at(i)) + std::to_string(assign_3.at(i));
			} else {
				path = std::to_string(assign_5.at(i)) + "-";
			}

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


// inspired by https://github.com/Eleobert/dbscan/blob/master/dbscan.cpp
std::vector<int> dbscan_aux(ClusterNode *curr_node, const int &points, const int &min_counts,
                            std::vector<std::vector<int>> &assignment, const bool &five) {


	int index;
	int clust_num = 0;
	std::vector<int> *adj_vec;
	std::vector<int> neighbors;
	std::vector<int> sub_neighbors;
	std::vector<int> assign_vec(points, -1);
	std::vector<bool> visted(points, false);

	// avoid copy
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
				if ((i != j) && (std::abs((*adj_vec)[j] - (*adj_vec)[i]) <= ImpaqtArguments::Args.epsilon)) {
					neighbors.push_back(j);
				}
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


void dbscan(ClusterList &cluster,  const int &strand) {

	int points;
	int min_counts;

	ClusterNode *curr_node = cluster.get_head(strand);

	while (curr_node != NULL) {

		points = curr_node -> get_read_count();

		if (points >= ImpaqtArguments::Args.min_count) {

			std::vector<std::string> paths;
			std::vector<int> assign_vec_5;
			std::vector<int> assign_vec_3;
			std::vector<std::vector<int>> assignments_5;
			std::vector<std::vector<int>> assignments_3;

			min_counts = std::max((int)((float)curr_node -> get_read_count() * ((float)ImpaqtArguments::Args.count_percentage / 100)), 20);

			// std::cerr << "///////////////////////////////////////////////\n";
			// std::cerr << "Region: " << curr_node -> get_chrom_index() << ":"
			//           << curr_node -> get_start() << "-"
			//           << curr_node -> get_stop() << "\n"
			//           << "Read Counts: " << curr_node -> get_read_count() << "\n"
			//           << "Min Count: " << min_counts << "\n"
			//           << "Epsilon: " << ImpaqtArguments::Args.epsilon
			//           << "\n///////////////////////////////////////////////\n";


			assign_vec_5 = dbscan_aux(curr_node, points, min_counts, assignments_5, true);
			assign_vec_3 = dbscan_aux(curr_node, points, min_counts, assignments_3, false);

			// If clusters were not found
			if (assignments_5.empty() && assignments_3.empty()) {
				curr_node = curr_node -> get_next();
				continue;
			}

			trace_transcripts(curr_node, paths, assign_vec_5, assign_vec_3);
			reduce_transcripts(curr_node, paths, assignments_5, assignments_3);

		}

		curr_node = curr_node -> get_next();
	}
}