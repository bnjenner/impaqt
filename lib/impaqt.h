#include "cluster.h"
#include "dbscan.h"
#include "annotation.h"
#include "api/BamAux.h"
#include "api/BamReader.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Impaqt Process Class
class Impaqt {

private:

	// Alignment file Readers
	BamTools::BamReader inFile;		   	 				 // Bam File Object
	BamTools::BamAlignment alignment;	  				 // BamAlignmentRecord record;

	static AnnotationList annotation;			  		 // Annotation
	static std::string alignment_file_name;	     		 // alignment file
	static std::string index;				  	     	 // alignment index file
	int chrom_index;					  				 // chromosome number
	bool ignore_chr = false;			  				 // to ignore for downstream

	static std::unordered_map<int, std::string> contig_map;
	static std::unordered_map<int, int> contig_lengths;

	size_t total_reads = 0;
	size_t unique_reads = 0;
	size_t ambiguous_reads = 0;
	size_t multimapped_reads = 0;
	size_t unassigned_reads = 0;


public:

	// Data Structures
	ClusterList cluster_list;			  // List for clusters

	/////////////////////////////////////////////////////////////
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

	/////////////////////////////////////////////////////////////
	// Get Reads Stats
	size_t get_total_reads() { return total_reads; }
	size_t get_unique_reads() { return unique_reads; }
    size_t get_ambiguous_reads() { return ambiguous_reads; }
    size_t get_multimapped_reads() { return multimapped_reads; }
    size_t get_unassigned_reads(){ return unassigned_reads; }

	// Get annotation
	AnnotationList* get_annotation() { return &annotation; }

	// Get Chromosome Info
	int get_chrom_num() { return contig_map.size(); }
	std::unordered_map<int, std::string> get_contig_map() { return contig_map; }
	std::unordered_map<int, int> get_contig_lengths() { return contig_lengths; }

	/////////////////////////////////////////////////////////////
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

	/////////////////////////////////////////////////////////////
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

	// Set Annotation (this might cause a memory leak)
	void add_annotation() {
		annotation = AnnotationList(); 
		annotation.create_gene_graph();
	}

	/////////////////////////////////////////////////////////////
	// Grab Alignments within Interval Using Bam Index
	void find_clusters() {

		cluster_list.initialize(chrom_index, contig_map[chrom_index], contig_lengths[chrom_index]);

		if (!inFile.Jump(chrom_index)) {
			std::cerr << "[ERROR: Could not jump to region: " << chrom_index << ".]\n";
			return;
		}

		// If failed to create clusters, flag to ignore
		if (!cluster_list.create_clusters(inFile, alignment)) { ignore_chr = true; }
	}

	// Merge neighboring clusters and remove zeroes
	void collapse_clusters() {
		cluster_list.collapse_clusters(0); // Forward
		cluster_list.collapse_clusters(1); // Reverse
}

	// Differentiate Transcripts
	void find_transcripts() {
		dbscan(cluster_list, 0); // Forward
		dbscan(cluster_list, 1); // Reverse
	}

	// Assign Transcripts to Genes
	void assign_transcripts() {
		// The work just actually has to be done here I think...
	}


	/////////////////////////////////////////////////////////////
	// Launch thread
	void launch() {
		this -> open_alignment_file();			  // open files
		this -> find_clusters();	  			  // find clusters
		if (ignore_chr) { return; }
		this -> collapse_clusters();	  		  // collapse clusters
		this -> find_transcripts();	  		  	  // dbscan clustering algorithm
		// this -> assign_transcripts();  	  			  // overlap genes
	}
	/////////////////////////////////////////////////////////////
};

