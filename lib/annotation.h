//////////////////////////////////////
// Annotation Class (doubley linked list)
class AnnotationList {

public:

	std::string feature_tag;
	std::string feature_id;
	std::string annotation_file_name;
	std::string chrom;
	std::string stranded;
	bool isGFF;

	GeneNode *head = NULL;
	GeneNode *tail = NULL;

	// Empty
	AnnotationList() {};

	// Proper constructor
	AnnotationList(const ImpaqtArguments *args) {
		annotation_file_name = args -> annotation_file;
		feature_tag = args -> feature_tag;
		feature_id = args -> feature_id;
		stranded = args -> stranded;
		isGFF = args -> isGFF;
	}

	// Copy constructor
	AnnotationList(const AnnotationList &anno) {
		annotation_file_name = anno.annotation_file_name;
		feature_tag = anno.feature_tag;
		head = anno.head;
		tail = anno.tail;
	}

	// Print genes and counts
	void print_genes() {
		GeneNode *curr_node = head;
		while (curr_node != NULL) {
			if (curr_node -> get_chrom() == chrom) {
				std::cout << curr_node -> get_gene() << "\t" << curr_node -> read_count << "\n";
			}
			curr_node = curr_node -> next;
		}
	}

	// Print annotation graph
	void print_graph() {
		GeneNode *curr_node = head;
		while (curr_node != NULL) {
			std::cout << curr_node -> gene_id << "\t";
			for (int i = 0; i < curr_node -> clust_count; i ++) {
				std::cout << curr_node -> clust_vec[(2 * i)] << "\t" << curr_node -> clust_vec[(2 * i) + 1] << "\t";
			}
			std::cout << "\n";
			curr_node = curr_node -> next;
		}
	}

	// Parse lines of annotation file
	void parse_annotation_line(const std::string &line, std::vector<std::string> &columns) {
		size_t n = 0;
		std::istringstream iss(line);
		std::string column;
		while (std::getline(iss, column, '\t')) { 
			columns.at(n) = column;
			n++;
		}
	}

	// Parse 9th column of GTF (tags)
	void parse_annotation_tags(const std::string &tag_column, std::vector<std::string> &tags, const bool &isGFF) {

		size_t n = 0;
		std::istringstream iss(tag_column);
		std::string tag;

		// GFF
		if (isGFF) {

			std::string subtag;

			while (std::getline(iss, tag, ';')) {
				std::istringstream iss2(tag);
				while (std::getline(iss2, subtag, '=')) { tags.push_back(subtag); }
			}

			// GTF
		} else {

			while (std::getline(iss, tag, ' ')) { 
				if (n % 2 == 1) {
					tag.resize(tag.size() - 2);
					tag.erase(0, 1);
				}
				tags.push_back(tag);
				n++;
			}
		}


	}

	// Create graph structure
	void create_gene_graph() {

		std::ifstream infile(annotation_file_name);

		if (!infile) {
			std::cerr << "ERROR: Could not read annotation file: " << annotation_file_name << "\n";
			throw "ERROR: Could not read annotation file.";
		}

		// Cceate variables for lines and annotation tags
		std::string line;
		std::vector<std::string> tags;
		std::vector<std::string> columns{9, ""};

		// Init temp variables
		int start;
		int stop;
		int temp_strand;
		std::string curr_id;
		std::string temp_id;
		
		// Initilize pointer
		GeneNode *curr_node = NULL;

		// Iterate through lines in file
		while (std::getline(infile, line)) {

			if (line.find_first_of('#') != 0) {

				parse_annotation_line(line, columns);

				if (columns.at(2).compare(feature_tag) == 0) {	// if correct type

					tags.clear();
					parse_annotation_tags(columns.at(8), tags, isGFF);

					// get id 
					for (int i = 0; i < tags.size(); i++) {
						if (tags[i] == feature_id) {
							temp_id = tags[i + 1];
							break;
						}
					}

					// if no feature ID field, throw error
					if (temp_id.empty()) {
						std::cerr << "ERROR: Could not identify Feature ID.\n";
						throw "ERROR: Could not identify Feature ID.";
					}

					start = std::stoi(columns.at(3)) - 1;
					stop = std::stoi(columns.at(4)) - 1;

					// if new feature
					if (temp_id.compare(curr_id) != 0) {

						curr_id = temp_id;

						if (stranded.compare("reverse") == 0) {
							temp_strand = 1 - ((columns.at(6) == "+") ? 0 : 1);
						} else {
							temp_strand = (columns.at(6) == "+") ? 0 : 1;
						}

						GeneNode *new_node = new GeneNode(temp_id,								// Gene ID
						    	                      	  temp_strand,							// Strand
						        	                 	  start, 								// Start
						            	        		  stop,   								// Stop
						                	       	  	  columns.at(0));						// Chrom

						// if first node
						if (curr_node == NULL) {
							curr_node = new_node;
							tail = curr_node;
							head = curr_node;

						} else {
							curr_node -> set_next(new_node);
							new_node -> set_prev(curr_node);
							curr_node = new_node;
							tail = curr_node;
						}


					} else {
						curr_node -> modify_cluster(std::stoi(columns.at(3)) - 1, 
													std::stoi(columns.at(4)) - 1, 
													0);
					}

				}
			}
		}

	}
};