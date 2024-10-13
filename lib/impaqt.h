#include "dbscan.h"

//////////////////////////////////////
// Impaqt Process Class
class Impaqt {

public:

	const ImpaqtArguments *parameters; 	  // parameters struct (found in parser.h)

	// Alignment file Readers
	BamTools::BamReader inFile;		   	  // Bam File Object
	BamTools::BamAlignment alignment;	  // BamAlignmentRecord record;

	// Data Structures
	ClusterList cluster_list;			  // List for clusters
	ClusterList transcript_list;		  // List for Transcripts 

	AnnotationList annotation;			  // Annotation 
	std::string alignment_file_name;	  // alignment file
	std::string index;				  	  // alignment index file
	int chr_num;						  // chromosome number
	bool ignore_chr = false;			  // to ignore for downstream 

	std::unordered_map<int, std::string> contig_map;

	size_t total_reads = 0;
	size_t ambiguous_reads = 0;
	size_t unique_reads = 0;
	size_t multimapped_reads = 0;
	size_t unassigned_reads = 0;

	// Empty
	Impaqt() {};

	// Initialized
	Impaqt(const ImpaqtArguments *args, int ref) {
		alignment_file_name = args -> alignment_file;
		index = args -> index_file;
		parameters = args;
		chr_num = ref;
	}


	void open_alignment_file() {
		// Open alignment file
		if (!inFile.Open(alignment_file_name)) {
			std::cerr << "ERROR: Could not read alignment file: " << alignment_file_name << "\n";
			throw "ERROR: Make sure alignment file exists.";
		}
		// Open index file
		if (!inFile.OpenIndex(index)) {
			std::cerr << "ERROR: Could not read index file: " << index << "\n";
			throw "ERROR: Make sure index is present in BAM file location.";
		}
	}

	void close_alignment_file() { inFile.Close(); }
	void print_counts() { annotation.print_genes(); }
	void print_gtf() { cluster_list.print_clusters(); }


	// parse input file for contig order and jump position
	void set_chrom_order() {

		BamTools::SamHeader head = inFile.GetHeader();

		if (head.HasSortOrder()) {
			std::string sortOrder = head.SortOrder;
			if (sortOrder.compare("coordinate") != 0) {
				std::cerr << "ERROR: Sorted alignment file required.\n";
				throw "ERROR: Could not read alignment file.";
			}
		} else {
			std::cerr << "ERROR: BAM file has no @HD SO:<SortOrder> attribute. Impossible to determine sort order.\n";
			throw "ERROR: Could determine sort status. Please ensure file is sorted.";

		}

		// Generate Ref Map (contig indicies)
		BamTools::RefVector references = inFile.GetReferenceData();
		for (int i = 0; i < references.size(); i++) {
			contig_map[i] = references.at(i).RefName;
		}
	}

	// Copy contig cache
	void copy_order(const std::unordered_map<int, std::string> &init_contig_map) {
		contig_map = init_contig_map;
	}

	// Copy annotation
	void copy_annotation(AnnotationList &init_annotation, const int &temp_chr_num) {
		/*
			Believe it or not, copying this is faster and the memory usage difference
				is negligible. I have no idea why this is the case.
		*/
		annotation = init_annotation;
		annotation.chrom = contig_map[temp_chr_num];
	}


	// Grab Alignments within Interval Using Bam Index
	void find_clusters() {

		cluster_list.initialize(chr_num, contig_map[chr_num], parameters);

		if (!inFile.Jump(chr_num)) {
			std::cerr << "[ERROR: Could not jump to region: " << chr_num << ".]\n";
			return;
		}

		if (cluster_list.set_head(inFile, alignment)) {
			// std::cerr << "head set!" << std::endl;
			cluster_list.create_clusters(inFile, alignment);  
			// std::cerr << "clusters found!" << std::endl;
			cluster_list.collapse_clusters();
			// std::cerr << "graph collaped!" << std::endl;

			multimapped_reads = cluster_list.multimapped_reads;
			total_reads = cluster_list.total_reads;
		
		} else {
			// std::cerr << "ignoring chr" << "\n";
			ignore_chr = true;
		}
	}

	void find_transcripts() {

		ClusterNode *curr_clust = cluster_list.head;
		ClusterNode *new_node;
		transcript_list.initialize(chr_num, contig_map[chr_num], parameters);

		while (curr_clust != NULL) {

			if (curr_clust -> read_count >= parameters -> min_count) {
				dbscan(curr_clust, transcript_list,
					   parameters -> count_percentage, 
					   parameters -> epsilon);
			} else {
				new_node = new ClusterNode(curr_clust -> clust_vec, curr_clust-> chrom_index, curr_clust -> strand, curr_clust -> read_count);
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

			// std::cerr << contig_map[chr_num] << " "
			// 	      << transcript_list.tail -> read_count << " "
			// 	      << transcript_list.tail -> clust_vec.at(0) << "\t";

			// for (int i = 0; i < ((transcript_list.tail -> clust_vec.size()) / 2); i++) {
			// 	std::cerr << transcript_list.tail -> clust_vec[(2*i)] << "-"
			// 			  << transcript_list.tail -> clust_vec[(2*i) + 1] << "\t";
			// }
			// std::cerr << "\n";

			curr_clust = curr_clust -> next;
		}

	}


	// Overlap genes
	void overlap_genes() {

		/*
			This function is horrible... I really need to
			find a better way to do this...
		*/

		if (cluster_list.head == NULL) { return; }

		// Graph pointers
		ClusterNode *curr_clust = cluster_list.head;
		GeneNode *curr_gene = annotation.head;
		GeneNode *temp_gene = NULL;

		// Vector of Pointers for start positions in list
		std::vector<GeneNode *> prev_gene{annotation.head, annotation.head};

		// Temp variables for overlap type and overlapping reads
		int overlap = 0;
		int temp_overlap = 0;
		int cmp_overlap = 0;
		int overlapping_reads = 0;
		int temp_overlapping_reads = 0;

		// Temp variables
		int temp_strand = 0;
		bool assigned = false;

		while (curr_clust != NULL) {

			assigned = false;

			// get strand, use it to index current gene to check
			temp_strand = curr_clust -> strand;
			curr_gene = prev_gene[temp_strand];


			// iterate through genes
			while (curr_gene != NULL) {

				// if on the wrong chromosome
				if (curr_gene -> chrom != contig_map[chr_num]) {

					// iterate through and find genes on correct contig
					while (curr_gene -> chrom != contig_map[chr_num])  {
						curr_gene = curr_gene -> next;
						if (curr_gene == NULL) {
							break;
						}
					}

					// iterate through and find genes on same strand
					while (curr_gene != NULL) {
						if (curr_gene -> strand == temp_strand) {
							break;
						}
						curr_gene = curr_gene -> next;
					}

					// Set gene starting point to current gene.
					prev_gene[temp_strand] = curr_gene;


					// gene is on correct contig
				} else {

					// If cluster is past gene, find next gene
					if (curr_clust -> get_start() > curr_gene -> get_stop()) {

						curr_gene = curr_gene -> next;

						// find next gene on same strand
						while (curr_gene != NULL) {
							if (curr_gene -> strand == temp_strand) {
								break;
							}
							curr_gene = curr_gene -> next;
						}


						// If cluster is before gene, find next cluster
					} else if (curr_clust -> get_stop() < curr_gene -> get_start()) {
						break;

						// if possibility of overlap
					} else {

						overlap = 0;
						overlapping_reads = 0;

						// iterate through subclusters. check overlap
						for (int x = 0; x < curr_clust -> clust_count; x++) {

							temp_overlap = curr_gene -> check_gene_overlap(curr_clust -> clust_vec[(2 * x)],
							                                        	   curr_clust -> clust_vec[(2 * x) + 1],
							                                        	   curr_clust -> strand);

							// if subcluster overlaps, add number of reads to total overlapping reads
							if (temp_overlap != 0) {
								overlapping_reads += curr_clust -> count_vec[x];
							}

							// check if max overlap
							overlap = std::max(temp_overlap, overlap);
						}


						// Check for ambigous overlaps (following genes)
						if (overlap != 0) {

							curr_clust -> assigned_gene	= curr_gene -> gene_id;
							temp_gene = curr_gene -> next;

							// Temp gene isn't NULL
							if (temp_gene != NULL) {

								// iterate through following genes
								while (temp_gene -> get_start() < curr_clust -> get_stop()) {

									// zero temp variables
									temp_overlapping_reads = 0;
									cmp_overlap = 0;

									// check for overlap, if there is, mark cluster as ambigous
									for (int x = 0; x < curr_clust -> clust_count; x++) {

										temp_overlap = temp_gene -> check_gene_overlap(curr_clust -> clust_vec[(2 * x)],
										                                        	   curr_clust -> clust_vec[(2 * x) + 1],
										                                        	   curr_clust -> strand);

										// if overlap, add to overlapping read accumulator
										if (temp_overlap > 0) {
											temp_overlapping_reads += curr_clust -> count_vec[x];
										}

										// keep track of max overlap
										cmp_overlap = std::max(temp_overlap, cmp_overlap);

									}

									// if number of overlapping reads is equal
									if (temp_overlapping_reads == overlapping_reads) {

										// if overlapping scores are equal, assign ambiguous
										if (cmp_overlap == overlap) {
											curr_clust -> ambiguous = 1;
											curr_clust -> assigned_gene	= "";
											break;

											// if new overlap is better, replace gene assignment and continue
										} else if (cmp_overlap > overlap) {
											curr_gene = temp_gene;
											overlapping_reads = temp_overlapping_reads;
											curr_clust -> assigned_gene	= curr_gene -> gene_id;
										}

										// if number of overlapping reads is greater than, replace gene assignment and continue
									} else if (temp_overlapping_reads > overlapping_reads) {
										curr_gene = temp_gene;
										overlapping_reads = temp_overlapping_reads;
										curr_clust -> assigned_gene	= curr_gene -> gene_id;
									}

									// move on to next gene
									temp_gene = temp_gene -> next;

									// if no more genes, break
									if (temp_gene == NULL) {
										break;
									}
								}
							}

							assigned = true;
							curr_clust -> unassigned = !assigned;

							// if not ambigous, assign reads to gene
							if (!(curr_clust -> ambiguous)) {
								curr_gene -> read_count += curr_clust -> read_count;
								curr_clust -> assigned_gene = curr_gene -> gene_id;
								unique_reads += curr_clust -> read_count;
							} else {
								ambiguous_reads += curr_clust -> read_count;
							}

							break;

						} else {

							// Issue is that gapped clusters do not advance cluster, yet skip many genes
							// in the process, solution is to trackback to genes

							curr_gene = curr_gene -> next;

							// find gene with correct strand
							while (curr_gene != NULL) {
								if (curr_gene -> strand == temp_strand) {
									break;
								}
								curr_gene = curr_gene -> next;
							}

						}

					}

				}

			}

			// increamenet unassigned reads
			if (assigned != true) {
				unassigned_reads += curr_clust -> read_count;
			}

			curr_clust = curr_clust -> next;
		}

	}

	void launch() {
		this -> open_alignment_file();			  // open files
		this -> find_clusters();	  			  // find clusters
		if (ignore_chr) { return; }
		this -> find_transcripts();	  		  // dbscan clustering algorithm
		this -> overlap_genes();  	  			  // overlap genes
		this -> close_alignment_file();		  	  // close files
	}
};

