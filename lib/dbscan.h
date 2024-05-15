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

void dbscan(ClusterNode *cluster, const int &count_percentage, const int &epsilon) {

	// inspired by https://github.com/Eleobert/dbscan/blob/master/dbscan.cpp

	/*
		IMPORTANT
			Creating a copy of the mid vector allows us to retrieve the index of the read mid point
			and therefore relate the assignment back to a read cluster, helping us establish read 
			cluster bounds.

			Should splices even be considered for clustering? If anything, it leadst to more
			ambiguity in their placement. Can clustering be considered a reliable quantification method?
			Ultimately my code would be 10x cleaner if I did not have to deal with them and they were
			instead considered individual reads. Like genuinely, they kinda just suck and can erroneously
			merge read clusters together when we don't want them to.
	*/

	int index;
	int min_counts = std::min((int)((float)cluster -> read_count * ((float)count_percentage / 100)), 20); 
	int points = cluster -> five_vec.size();
	std::vector<bool> visted(points, false);
	std::vector<std::vector<int>> assignment;
	std::vector<int> neighbors;
	std::vector<int> sub_neighbors;

	std::vector<int> adj_five_vec = offset(cluster);

	std::cerr << "////////////////////////////////////////////////\n";
	std::cerr << "Region: " << cluster -> chrom_index << ":"
			  << cluster -> get_start() << "-"
			  << cluster -> get_stop() << "\n"
			  << "Read Counts: " << cluster -> read_count << "\n"
			  << "Min Count: " << min_counts << "\n"
			  << "Epsilon: " << epsilon; 

	// for (const auto &clust : cluster -> clust_vec) {
	// 	std::cerr << clust << "\t";
	// }

	std::cerr << "\n///////////////////////\n";

	for (const auto &mid : adj_five_vec) {
		std::cerr << mid << "\t";
	}

	std::cerr << "\n///////////////////////\n";

	for (int i = 0; i < points; i++) {

		neighbors.clear();

		if (visted[i] == false) {

			for (int j = 0; j < points; j++) {
				if ((i != j) && (std::abs(adj_five_vec[j] - adj_five_vec[i]) <= epsilon)) {
					neighbors.push_back(j);
					// std::cerr << adj_five_vec[i] << "\t" << adj_five_vec[j] << "\n\t"
					// 		  << std::abs(adj_five_vec[i] - adj_five_vec[j]) << "\n";
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


	int min_offset;
	int max_offset;
	std::vector<int>::iterator min_result;
	std::vector<int>::iterator max_result;

	for (int i = 0; i < assignment.size(); i++) {

		if (cluster -> strand == 0) {
			min_offset = -epsilon;
			max_offset = 0;
		} else {
			min_offset = 0;
			max_offset = epsilon;
		}

		std::cerr << "Cluster: " << i << "\n";
		std::cerr << "    Core Points: " << assignment[i].size() << "\n\t";
		min_result = std::min_element(assignment[i].begin(), assignment[i].end());
		max_result = std::max_element(assignment[i].begin(), assignment[i].end());
		std::cerr << cluster -> five_vec.at(*min_result) + min_offset << "-" 
				  << cluster -> five_vec.at(*max_result) + max_offset << "\n\t";
		std::cerr << adj_five_vec.at(*min_result) << "-" 
				  << adj_five_vec.at(*max_result) << "\n";

		std::cerr << "///////////////////////\n";

	}

}