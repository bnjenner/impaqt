//////////////////////////////////////
// Annotation Class (doubley linked list)
class AnnotationList {

public:

	std::string annotation_file_name;
	std::string chrom;
	std::string stranded;
	bool isGFF;

	std::string feature_id;
	std::string feature_tag;

	// GeneNode *head = NULL;
	// GeneNode *tail = NULL;

	// Empty
	AnnotationList() {};

	// Proper constructor
	AnnotationList(const ImpaqtArguments *args) {
		annotation_file_name = args -> annotation_file;
		feature_id = args -> feature_id;
		feature_tag = args -> feature_tag;
		stranded = args -> stranded;
		isGFF = args -> isGFF;
	}

	// Copy constructor
	AnnotationList(const AnnotationList &anno) {
		annotation_file_name = anno.annotation_file_name;
		feature_id = anno.feature_id;
		feature_tag = anno.feature_tag;
		stranded = anno.stranded;
		isGFF = anno.isGFF;
		// head = anno.head;
		// tail = anno.tail;
	}

	// // Print genes and counts
	// void print_genes() {
	// 	GeneNode *curr_node = head;
	// 	while (curr_node != NULL) {
	// 		if (curr_node -> get_chrom() == chrom) {
	// 			std::cout << curr_node -> get_gene() << "\t" << curr_node -> read_count << "\n";
	// 		}
	// 		curr_node = curr_node -> next;
	// 	}
	// }

	// // Parse lines of annotation file
	// void parse_annotation_line(const std::string &line, std::vector<std::string> &columns) {
	// 	size_t i = 0;
	// 	std::istringstream iss(line);
	// 	std::string column;
	// 	while (std::getline(iss, column, '\t')) {
	// 		columns.at(i) = column;
	// 		i++;
	// 	}
	// }

	// // Get feature ID
	// std::string get_feature_id(const std::string &anno_column, const bool &isGFF) {

	// 	std::string tag;
	// 	std::string subtag;

	// 	char sep = '"';
	// 	size_t i = 0;

	// 	std::istringstream iss(anno_column);

	// 	if (isGFF) { sep = '='; }

	// 	while (std::getline(iss, tag, ';')) {

	// 		if (tag.find(feature_id) != std::string::npos) {
	// 			std::istringstream iss2(tag);

	// 			while (std::getline(iss2, subtag, sep)) {
	// 				if (i % 2 == 1) {
	// 					return subtag;
	// 				}
	// 				i++;
	// 			}
	// 		}
	// 	}

	// 	throw "ERROR: Could not find feature tag in line of annotation file. Check consistency of formatting.";
	// }

	// // Create graph structure
	// void create_gene_graph() {

	// 	std::ifstream infile(annotation_file_name);

		// if (!infile) { throw "ERROR: Could not read annotation file."; }

		// // Cceate variables for gtf columns
		// std::string line;
		// std::vector<std::string> columns{9, ""};

		// // Init temp variables
		// int start;
		// int stop;
		// int temp_strand;
		// std::string curr_id;
		// std::string temp_id;

		// // Initilize pointer
		// GeneNode *curr_node = NULL;

		// // Iterate through lines in file
		// while (std::getline(infile, line)) {

		// 	if (line.find_first_of('#') != 0) {
		// 		parse_annotation_line(line, columns);

		// 		if (columns.at(2).compare(feature_tag) == 0) {
		// 			temp_id = get_feature_id(columns.at(8), isGFF);

		// 			// if no feature ID field, throw error
	// 				if (temp_id.empty()) {
	// 					std::cerr << "ERROR: Could not identify Feature ID.\n";
	// 					throw "ERROR: Could not identify Feature ID.";
	// 				}

	// 				start = std::stoi(columns.at(3)) - 1;
	// 				stop = std::stoi(columns.at(4)) - 1;

	// 				// if new feature
	// 				if (temp_id.compare(curr_id) != 0) {

	// 					curr_id = temp_id;

	// 					if (stranded.compare("reverse") == 0) {
	// 						temp_strand = 1 - ((columns.at(6) == "+") ? 0 : 1);
	// 					} else {
	// 						temp_strand = (columns.at(6) == "+") ? 0 : 1;
	// 					}

	// 					GeneNode *new_node = new GeneNode(temp_id,			// Gene ID
	// 					                                  temp_strand,		// Strand
	// 					                                  start, 			// Start
	// 					                                  stop,   			// Stop
	// 					                                  columns.at(0));	// Chrom

	// 					// if first node
	// 					if (curr_node == NULL) {
	// 						curr_node = new_node;
	// 						tail = curr_node;
	// 						head = curr_node;

	// 					} else {
	// 						curr_node -> set_next(new_node);
	// 						new_node -> set_prev(curr_node);
	// 						curr_node = new_node;
	// 						tail = curr_node;
	// 					}


	// 				} else {
	// 					curr_node -> modify_cluster(std::stoi(columns.at(3)) - 1,
	// 					                            std::stoi(columns.at(4)) - 1,
	// 					                            0);
	// 				}

	// 			}
	// 		}
	// 	}

	// }
};