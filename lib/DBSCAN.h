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


// Merge Overlapping Transcripts
void reduce_transcripts(ClusterNode *curr_node, std::vector<std::vector<int>> &transcripts, std::vector<float> &counts) {

	int pos;

	// If reverse strand, reverse and make negative (I'm actually pretty proud of this solution)
	if (curr_node -> get_strand() == 1) {
		for (int i = 0; i < transcripts.size(); i++) {
			transcripts[i] = reverse_and_negate(transcripts[i]);
		}
		std::reverse(transcripts.begin(), transcripts.end());
	}

	for (int i = 0; i < transcripts.size(); i++) {

		for (int j = i + 1; j < transcripts.size(); j++) {
			pos = transcripts[i].back();

			if (pos >= transcripts[j][0] && pos <= transcripts[j][1]) {
				transcripts[i].back() = transcripts[j][1];

				if (transcripts[j].size() > 2) {
					transcripts[i].push_back(transcripts[j][2]);
					transcripts[i].push_back(transcripts[j][3]);
				}

				transcripts[j].clear();
				i = j;
			}

			j += 1;
		}

		if (!transcripts[i].empty()) {

			// Reinstate original order
			if (curr_node -> get_strand() == 1) {
				transcripts[i] = reverse_and_negate(transcripts[i]);
			}

			curr_node -> add_transcript(transcripts[i], counts[i]);
		}
	}
}


// Get Transcript Coordinates
void get_transcripts(ClusterNode *curr_node, const std::vector<std::string> &paths,
                     	std::vector<std::vector<int>> &core_5, std::vector<std::vector<int>> &core_3,
                     	std::vector<std::vector<int>> *transcripts, std::vector<float> *counts) {

	float core_points;
	int index, clusters;
	int min_pos, max_pos;  // necessary because 3' vec may not be in order due to splicing

	// Get transcript coordinates
	//	BNJ: 5/2/2025 - Horrible extending vector, but they're short so it shouldn't matter THAT much
	//  BNJ: 5/30/202 - Also the damn minus strand man. Idk why it's never in order. 

	for (const auto &p : paths) {

		std::cerr << "\t" << p << "\n";

		core_points = 0.0f;
		clusters = 0;
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
			clusters += 1;

			std::cerr << "\t" << tmp_vec[0] << "-" << tmp_vec[1] << "\n";
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
			clusters += 1;

			std::cerr << "\t" << min_pos << "-" << max_pos << "\n";
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
						   	   tmp_vec[0], tmp_vec[1]};

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
		std::cerr << "\t" << p << " : ";
		for (const auto &t : tmp_vec) { std::cerr << t << " "; }
		std::cerr << "\n";
		counts -> push_back(core_points / clusters);
	}
}


// Find all linked DBSCAN clusters
void trace_transcripts(ClusterNode *curr_node, std::vector<std::string> &path_vec,
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

	/*

	OK SO, the merging of paths needs to happen prior to converting the indexes to coordinates. thats the ticket. 

	*/



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
			// min_counts = (int)((float)points * ((float)ImpaqtArguments::Args.count_percentage / 100));

			std::cerr << curr_node -> get_contig_name() << ":" 
			          << curr_node -> get_start() << "-" << curr_node -> get_stop() 
			          << " (" << points << " / " << min_counts << "\n";

			// Run DBSCAN
			assign_vec_5 = dbscan(curr_node, points, min_counts, assignments_5, true);
			assign_vec_3 = dbscan(curr_node, points, min_counts, assignments_3, false);
			

			// If clusters were not found
			if (assignments_5.empty() && assignments_3.empty()) {
				curr_node = curr_node -> get_next();
				continue;
			}

			// Find all linked DBSCAN clusters 
			trace_transcripts(curr_node, paths, assign_vec_5, assign_vec_3);	  		 
			
			// Get Transcript Coords and Core Points
			get_transcripts(curr_node, paths, 
							assignments_5, assignments_3, 
							&transcripts, &counts); 
			
			// Merge overlapping transcripts and add to cluster node
			reduce_transcripts(curr_node, transcripts, counts);	
		}

		curr_node = curr_node -> get_next();
	}
}