// inspired by https://github.com/Eleobert/dbscan/blob/master/dbscan.cpp
void dbscan_aux(ClusterNode *curr_node, const int &points, const int &min_counts, const int &epsilon, bool five) {

	/*
	So here is what we are going to do:
		Cluster by 5' ends and then cluster by 3' ends
		See what clusters overlap within an readlength between 5 and 3 clusters
		
			- If 1 cluster is identified and they overlap within the read length, they are one simple clust
			- If N clusters are identified each, and their overlap with their other end counterpart, we have N simple clusters
	*/

	int index;
	std::vector<bool> visted(points, false);
	std::vector<std::vector<int>> assignment;
	std::vector<int> neighbors;
	std::vector<int> sub_neighbors;

	std::vector<int>::iterator min_result;
	std::vector<int>::iterator max_result;	

	std::vector<int> adj_vec;

	if (five) {
		adj_vec = curr_node -> get_five_vec();
	} else {
		adj_vec = curr_node -> get_three_vec();
	}


	for (int i = 0; i < points; i++) {

		neighbors.clear();

		if (visted[i] == false) {

			for (int j = 0; j < points; j++) {
				if ((i != j) && (std::abs(adj_vec[j] - adj_vec[i]) <= epsilon)) {
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
							if (std::abs(adj_vec[index] - adj_vec[k]) <= epsilon) {
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



	std::cout << "/////////////////////\n";
	std::cout << "Chrom: " << curr_node -> get_chrom_index() << "\n";
	std::cout << "Five: " << five << "\n";
	std::cout << "Clusters Identified: " << assignment.size() << "\n";
	std::cout << "Regions: " << curr_node -> get_start() << "-" << curr_node -> get_stop() << "\n";
	std::cout << "//////////////////////\n";


	for (int i = 0; i < assignment.size(); i++) {

		std::cout << "Cluster: " << i << "\n";
		std::cout << "    Core Points: " << assignment[i].size() << "\n\t";
		min_result = std::min_element(assignment[i].begin(), assignment[i].end());
		max_result = std::max_element(assignment[i].begin(), assignment[i].end());
		std::cout << adj_vec.at(*min_result) << "-"
		          << adj_vec.at(*max_result) << "\n";

		std::cout << "///////////////////////\n";
	}

}


void dbscan(ClusterList &cluster,  const int &strand, const int &count_percentage, 
			const int &epsilon,  const int &min_count) {

	int points;
	int min_counts;

	ClusterNode *curr_node = cluster.get_head(strand);

	while (curr_node != NULL) {

		if (curr_node -> get_read_count() >= min_count) {

			points = curr_node -> get_read_count();
			min_counts = std::max((int)((float)curr_node -> get_read_count() * ((float)count_percentage / 100)), 20);

			std::cout << "////////////////////////////////////////////////\n";
			std::cout << "Region: " << curr_node -> get_chrom_index() << ":"
					  << curr_node -> get_start() << "-"
					  << curr_node -> get_stop() << "\n"
					  << "Read Counts: " << curr_node -> get_read_count() << "\n"
					  << "Min Count: " << min_counts << "\n"
					  << "Epsilon: " << epsilon 
					  << "\n///////////////////////////////////////////////\n";


			dbscan_aux(curr_node, points, min_counts, epsilon, true);
			dbscan_aux(curr_node, points, min_counts, epsilon, false);

		}

		curr_node = curr_node -> get_next();
	}
}