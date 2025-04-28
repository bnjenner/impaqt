#include "cluster.h"
#include "dbscan.h"
#include "api/BamAux.h"
#include "api/BamReader.h"

//////////////////////////////////////
// Impaqt Process Class
class Impaqt {

private:

	// Alignment file Readers
	BamTools::BamReader inFile;		   	 				 // Bam File Object
	BamTools::BamAlignment alignment;	  				 // BamAlignmentRecord record;

	// AnnotationList annotation;			  			 // Annotation
	static std::string alignment_file_name;	     		 // alignment file
	static std::string index;				  	     	 // alignment index file
	int chrom_index;					  				 // chromosome number
	bool ignore_chr = false;			  				 // to ignore for downstream

	static std::unordered_map<int, std::string> contig_map;
	static std::unordered_map<int, int> contig_lengths;

	size_t total_reads = 0;
	size_t ambiguous_reads = 0;
	size_t unique_reads = 0;
	size_t multimapped_reads = 0;
	size_t unassigned_reads = 0;


public:

	// Data Structures
	ClusterList cluster_list;			  // List for clusters
	// ClusterList transcript_list;		  // List for Transcripts

	// Empty
	Impaqt() {};

	// Initialized
	Impaqt(int ref) {
		alignment_file_name = ImpaqtArguments::Args.alignment_file;
		index = ImpaqtArguments::Args.index_file;
		chrom_index = ref;
	}

	// Empty
	~Impaqt() {
		close_alignment_file();
	};


	// Get Reads Stats
	size_t get_total_reads() { return total_reads; }
	size_t get_multimapped_reads() { return multimapped_reads; }
	std::unordered_map<int, std::string> get_contig_map() { return contig_map; }
	int get_chrom_num() { return contig_map.size(); }
	std::unordered_map<int, int> get_contig_lengths() { return contig_lengths; }

	// Open BAM file
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

	// Close Bam File
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

		cluster_list.initialize(chrom_index, contig_map[chrom_index], contig_lengths[chrom_index]);

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
		cluster_list.collapse_clusters(0); // Forward
		cluster_list.collapse_clusters(1); // Reverse
}

	// Differentiate Transcriptes
	void find_transcripts() {
		dbscan(cluster_list, 0);
		dbscan(cluster_list, 1);
	}


	void launch() {
		this -> open_alignment_file();			  // open files
		this -> find_clusters();	  			  // find clusters
		if (ignore_chr) { return; }
		this -> collapse_clusters();	  		  // collapse clusters
		this -> find_transcripts();	  		  	  // dbscan clustering algorithm
		// this -> overlap_genes();  	  			  // overlap genes
	}
};

