// inspired by https://github.com/Eleobert/dbscan/blob/master/dbscan.cpp
void dbscan_aux(ClusterNode *curr_node, const int &count_percentage, const int &epsilon) {

	int min_counts = std::min((int)((float)curr_node -> get_read_count() * ((float)count_percentage / 100)), 20);
	int points = curr_node -> get_read_count();

	int index;
	std::vector<bool> visted(points, false);
	std::vector<std::vector<int>> assignment;
	std::vector<int> neighbors;
	std::vector<int> sub_neighbors;

	std::vector<int>::iterator min_result;
	std::vector<int>::iterator max_result;	

	std::cerr << "////////////////////////////////////////////////\n";
	std::cerr << "Region: " << curr_node -> get_chrom_index() << ":"
			  << curr_node -> get_start() << "-"
			  << curr_node -> get_stop() << "\n"
			  << "Read Counts: " << curr_node -> get_read_count() << "\n"
			  << "Min Count: " << min_counts << "\n"
			  << "Epsilon: " << epsilon 
			  << "\n///////////////////////////////////////////////\n";

	std::vector<int> adj_five_vec = curr_node -> get_five_vec();

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



	std::cerr << "//////////////////////////////////////////////\n";
	std::cerr << "Chrom: " << curr_node -> get_chrom_index() << "\n";
	std::cerr << "Clusters Identified: " << assignment.size() << "\n";
	std::cerr << "Regions: " << curr_node -> get_start() << "-" << curr_node -> get_stop() << "\n";
	std::cerr << "///////////////////////\n";


	for (int i = 0; i < assignment.size(); i++) {

		std::cerr << "Cluster: " << i << "\n";
		std::cerr << "    Core Points: " << assignment[i].size() << "\n\t";
		min_result = std::min_element(assignment[i].begin(), assignment[i].end());
		max_result = std::max_element(assignment[i].begin(), assignment[i].end());
		std::cerr << curr_node -> get_five_vec().at(*min_result) << "-"
		          << curr_node -> get_five_vec().at(*max_result) << "\n";

		std::cerr << "///////////////////////\n";
	}

}


void dbscan(ClusterList &cluster,  const int &strand, const int &count_percentage, 
			const int &epsilon,  const int &min_count) {

	ClusterNode *curr_node = cluster.get_head(strand);

	while (curr_node != NULL) {

		if (curr_node -> get_read_count() >= min_count) {
			dbscan_aux(curr_node, count_percentage, epsilon);
		}
		
		curr_node = curr_node -> get_next();
	}
}