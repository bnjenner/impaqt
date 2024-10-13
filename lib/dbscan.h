std::vector<int> offset(ClusterNode *cluster) {

	std::vector<int> adj_five_vec = cluster -> five_vec;

	for (int i = 0; i < adj_five_vec.size(); i++) {
			for (int j = 1; j < cluster -> clust_count; j++) {
				if (adj_five_vec[i] >= cluster -> clust_vec[(2 * j)]) {
					adj_five_vec[i] -= (cluster -> clust_vec[(2 * j)] - cluster -> clust_vec[(2 * j) - 1]);
				} else {
					break;
				}
			}
			adj_five_vec[i] -= cluster -> clust_vec[0];
		}
		std::sort(adj_five_vec.begin(), adj_five_vec.end()); 

	return adj_five_vec;
}

void dbscan(ClusterNode *cluster, ClusterList &transcript_list, const int &count_percentage, const int &epsilon) {

	// inspired by https://github.com/Eleobert/dbscan/blob/master/dbscan.cpp


	int index;
	int min_counts = std::min((int)((float)cluster -> read_count * ((float)count_percentage / 100)), 20); 
	int points = cluster -> five_vec.size();
	std::vector<bool> visted(points, false);
	std::vector<std::vector<int>> assignment;
	std::vector<int> neighbors;
	std::vector<int> sub_neighbors;

	std::vector<int> adj_five_vec = offset(cluster);

	// std::cerr << "////////////////////////////////////////////////\n";
	// std::cerr << "Region: " << cluster -> chrom_index << ":"
	// 		  << cluster -> get_start() << "-"
	// 		  << cluster -> get_stop() << "\n"
	// 		  << "Cluster Counts: " << cluster -> count_vec.size() << "\n"
	// 		  << "Read Counts: " << cluster -> read_count << "\n"
	// 		  << "Min Count: " << min_counts << "\n"
	// 		  << "Epsilon: " << epsilon; 


	// std::cerr << "\n///////////////////////\n";



	// std::cerr << "\n///////////////////////\n";

	for (int i = 0; i < points; i++) {

		neighbors.clear();

		if (visted[i] == false) {

			for (int j = 0; j < points; j++) {
				if ((i != j) && (std::abs(adj_five_vec[j] - adj_five_vec[i]) <= epsilon)) {
					neighbors.push_back(j);
				}
			}

			if (neighbors.size() >= min_counts) {

				std::vector<int> cluster_indexes;
				visted[i] = true;

				while (neighbors.empty() == false) {

					sub_neighbors.clear();
					index = neighbors.back();
					neighbors.pop_back();

					if (visted[index] == false) {

						visted[index] = true;

						for (int k = 0; k < points; k++) {
							if (std::abs(adj_five_vec[index] - adj_five_vec[k]) <= epsilon) {
								sub_neighbors.push_back(k);
							}
						}

						if (sub_neighbors.size() >= min_counts) {
							std::copy(sub_neighbors.begin(), sub_neighbors.end(), std::back_inserter(neighbors));
						}

						cluster_indexes.push_back(index);

					}
				}

				assignment.emplace_back(std::move(cluster_indexes));
			}
		}
	}


	int min_offset = 0;
	int max_offset = 0;
	int temp_count = 0;
	std::vector<int> temp_vec;
	std::vector<int>::iterator min_result;
	std::vector<int>::iterator max_result;
	ClusterNode *new_node;

	/*
	So here is the deal (and it is a shitty deal)
		Assigning clusters to genes is gonna require some decisions on lots of possible scenarios
		To make this consistent, let's use some guiding principles. 

		1. ASSIGNMENT TO GENES SHOULD ONLY BE DETERMINED BY THE CLUSTERS
		2. CLUSTERS, IF SPANNING JUNCTIONS, SHOULD BE REPRESENTED AS SUBCLUSTERS
		3

		FUCK YOU I AM MAKING AN EXECUTIVE DECISION, IF I CAN'T RELIABLE ASSIGN IT TO A CLUSTER IT DOESNT FUCKING EXIST


	*/

	// std::cerr << "\t" << cluster -> chrom_index 
	// 		  << "\t" << cluster -> clust_vec[0]
	// 		  << "\t" << cluster -> clust_vec[cluster -> clust_vec.size() - 1]
			  // << "\n";

	for (int i = 0; i < assignment.size(); i++) {

		// std::cerr << "Cluster: " << i << "\n";
		// std::cerr << "    Core Points: " << assignment[i].size() << "\n\t";
		// min_result = std::min_element(assignment[i].begin(), assignment[i].end());
		// max_result = std::max_element(assignment[i].begin(), assignment[i].end());
		// std::cerr << cluster -> five_vec.at(*min_result) + min_offset << "-" 
		// 		  << cluster -> five_vec.at(*max_result) + max_offset << "\n\t";
		// std::cerr << adj_five_vec.at(*min_result) << "-" 
		// 		  << adj_five_vec.at(*max_result) << "\n";

		// std::cerr << "///////////////////////\n";

		min_result = std::min_element(assignment[i].begin(), assignment[i].end());
		max_result = std::max_element(assignment[i].begin(), assignment[i].end());

		temp_vec = {cluster -> five_vec.at(*min_result)};

		for (int j = 1; j < (cluster -> clust_vec.size()) - 1; j += 2) {

			if (cluster -> five_vec.at(*min_result) <= cluster -> clust_vec.at(j)) {

				// If transcript zone spans splice junction 
				if ((cluster -> five_vec.at(*max_result) > cluster -> clust_vec.at(j + 1))) {
					temp_vec.push_back(cluster -> clust_vec.at(j));
					temp_vec.push_back(cluster -> clust_vec.at(j + 1));
				
				} else if (cluster -> five_vec.at(*max_result) <= cluster -> clust_vec.at(j + 1)){
					temp_vec.push_back(cluster -> five_vec.at(*max_result));
					break;
				} else if (cluster -> five_vec.at(*max_result) <= cluster -> clust_vec.at(j)) {
					temp_vec.push_back(cluster -> five_vec.at(*max_result));
					break;
				}

			} 
		}

		if (temp_vec.size() == 1) {
			temp_vec.push_back(cluster -> five_vec.at(*max_result));
			break;
		}

		// std::vector<int> &temp_vec, int ref_num, int temp_strand, int temp_count
		new_node = new ClusterNode(temp_vec, cluster -> chrom_index, cluster -> strand, assignment[i].size());

		if (transcript_list.hashead == false) {
			transcript_list.head = new_node;
			new_node -> ishead = true;
			transcript_list.hashead = true;
			transcript_list.tail = new_node;
		
		} else {
			new_node -> set_prev(transcript_list.tail);
			transcript_list.tail -> set_next(new_node);
			transcript_list.tail = new_node;
		}

	}

}