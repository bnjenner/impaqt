#include "cluster.h"
#include "dbscan.h"
#include "api/BamAux.h"
#include "api/BamReader.h"

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
	// ClusterList transcript_list;		  // List for Transcripts

	AnnotationList annotation;			  // Annotation
	std::string alignment_file_name;	  // alignment file
	std::string index;				  	  // alignment index file
	int chrom_index;					  // chromosome number
	bool ignore_chr = false;			  // to ignore for downstream

	std::unordered_map<int, std::string> contig_map;
	std::unordered_map<int, int> contig_lengths;

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
		chrom_index = ref;
	}

	// Empty
	~Impaqt() {
		contig_map.clear();
		contig_lengths.clear();
		close_alignment_file();
	};


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
	// void print_counts() { annotation.print_genes(); }
	// void print_gtf() { if (!ignore_chr) { cluster_list.print_clusters(); } }

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
			contig_lengths[i] = references.at(i).RefLength;
		}
	}

	// Copy contig cache
	void copy_order(const std::unordered_map<int, std::string> &t_contig_map,
					const std::unordered_map<int, int> &t_contig_lengths) {
		// potentially needless copy, will address
		contig_map = t_contig_map;
		contig_lengths = t_contig_lengths;
	}

	// // Copy annotation
	// void copy_annotation(AnnotationList &t_annotation, const int &t_chrom_index) {
	// 	/*
	// 		Believe it or not, copying this is faster and the memory usage difference
	// 			is negligible. I have no idea why this is the case.

	// 			TO DO: Double check this
	// 	*/
	// 	annotation = t_annotation;
	// 	annotation.chrom = contig_map[t_chrom_index];
	// }


	// Grab Alignments within Interval Using Bam Index
	void find_clusters() {

		cluster_list.initialize(chrom_index, contig_map[chrom_index], contig_lengths[chrom_index], parameters);

		if (!inFile.Jump(chrom_index)) {
			std::cerr << "[ERROR: Could not jump to region: " << chrom_index << ".]\n";
			return;
		}

		if (cluster_list.create_clusters(inFile, alignment)) {
			multimapped_reads = cluster_list.get_multimapped_reads();
			total_reads = cluster_list.get_total_reads();

		} else {
			ignore_chr = true;
		}
	}

	// Merge neighboring clusters and remove zeroes
	void collapse_clusters() {
		if (!ignore_chr) {
			cluster_list.collapse_clusters(0); // Forward
			cluster_list.collapse_clusters(1); // Reverse		
		}
	}

	// Differentiate Transcriptes
	void find_transcripts() {
		if (!ignore_chr) {
			dbscan(cluster_list, 0,
		       parameters -> count_percentage,
		       parameters -> epsilon);
			dbscan(cluster_list, 1,
		       parameters -> count_percentage,
		       parameters -> epsilon);
		}
	}


	void launch() {
		this -> open_alignment_file();			  // open files
		this -> find_clusters();	  			  // find clusters
		this -> collapse_clusters();	  		  // collapse clusters
		this -> find_transcripts();	  		  	  // dbscan clustering algorithm
		// this -> overlap_genes();  	  			  // overlap genes
	}
};

