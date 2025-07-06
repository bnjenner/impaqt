#include <unordered_map>
#include <api/BamAux.h>
#include <api/BamReader.h>

#include "AnnotationList.h"
#include "ClusterList.h"
#include "DBSCAN.h"
#include "AssignClusters.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Impaqt Process Class */

class Impaqt {

private:

	// Alignment file Readers
	BamTools::BamReader inFile;                                 // Bam File Object
	BamTools::BamAlignment alignment;                           // BamAlignmentRecord Object;

	static AnnotationList annotation;                           // Annotation list for genes
	static std::string alignment_file_name;                     // alignment file
	static std::string index_file_name;                         // alignment index file

	// These could just be arrays or vectors
	static std::unordered_map<int, std::string> contig_map;     // Links Index to Contig Name
	static std::unordered_map<int, int> contig_lengths;         // Links Index to Contig Length
	
	int contig_index;                                           // chromosome number
	bool ignore = false;                                        // to ignore for downstream

	std::shared_ptr<ClusterList> list_ptr;                      // Temporary Pointer to ClusterList

	// Read Stats
	long double assigned_reads = 0.0;
	long double unassigned_reads = 0.0;
	long double ambiguous_reads = 0.0;
	size_t multimapped_reads = 0;
	size_t low_quality_reads = 0;
	size_t total_reads = 0;
	size_t transcript_num = 0;


public:

	/////////////////////////////////////////////////////////////
	/* Constructors */

	// Empty
	Impaqt() {};

	// Initialized
	Impaqt(const int &contig_index) {
		this -> contig_index = contig_index;
		this -> index_file_name = ImpaqtArguments::Args.index_file;
		this -> alignment_file_name = ImpaqtArguments::Args.alignment_file;
	}

	/////////////////////////////////////////////////////////////
	/* Get Functions */

	// Get Reads Stats
	long double get_assigned_reads() { return assigned_reads; }
	long double get_unassigned_reads() { return unassigned_reads; }
	long double get_ambiguous_reads() { return ambiguous_reads; }

	size_t get_multimapped_reads() { return multimapped_reads; }
	size_t get_low_quality_reads() { return low_quality_reads; }
	size_t get_total_reads() { return total_reads; }
	size_t get_transcript_num() { return transcript_num; }

	// Get Data Structures
	AnnotationList* get_annotation() { return &annotation; }

	// Get Chromosome Info
	bool is_ignored() { return ignore; }
	int get_chrom_index() { return contig_index; }
	int get_chrom_num() { return contig_map.size(); }
	int get_contig_length() { return contig_lengths[contig_index]; }
	std::string get_contig_name() { return contig_map[contig_index]; }	
	std::unordered_map<int, std::string> get_contig_map() { return contig_map; }
	std::unordered_map<int, int> get_contig_lengths() { return contig_lengths; }


	/////////////////////////////////////////////////////////////
	/* Thread Initilizers */

	void open_alignment_file() {
		// Open alignment file
		if (!inFile.Open(alignment_file_name)) {
			std::cerr << "ERROR: Could not read alignment file: " << alignment_file_name << "\n";
			throw "ERROR: Make sure alignment file exists.";
		}
		// Open index file
		if (!inFile.OpenIndex(index_file_name)) {
			std::cerr << "ERROR: Could not read index file: " << index_file_name << "\n";
			throw "ERROR: Make sure index is present in BAM file location.";
		}
	}

	void close_alignment_file() { inFile.Close(); }

	// Parse input file for contig order and jump positions
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

	// Add and Create GTF/GFF Gene Annotation
	void add_annotation() {
		annotation = AnnotationList();
		annotation.create_gene_list();
	}

	/////////////////////////////////////////////////////////////
	/* Cluster Related Functions */

	// Create cluster list
	void create_clusters(std::shared_ptr<ClusterList> &cluster) {

		cluster = std::make_shared<ClusterList>();
		cluster -> initialize(contig_index, this -> get_contig_name(), this -> get_contig_length());

		if (!inFile.Jump(contig_index)) {
			std::cerr << "//ERROR: Could not jump to region: " << this -> get_contig_name() << "\n";
			throw "ERROR: Could not jump to region. Make sure BAM header is correct.";
		}

		// If failed to create clusters, flag to ignore
		if (!(cluster -> create_clusters(inFile, alignment))) { ignore = true; }
	}

	// Merge neighboring clusters and remove zeroes
	void collapse_clusters(std::shared_ptr<ClusterList> &cluster) {
		cluster -> collapse_clusters(0); // Forward 
		cluster -> collapse_clusters(1); // Reverse
	}

	// Differentiate Transcripts
	void find_transcripts(std::shared_ptr<ClusterList> &cluster) {
		if (ignore) { return; }

		// Forward
		find_transcripts_DBSCAN(cluster, 0);
		collapse_final_transcripts(cluster, 0);

		// Reverse
		find_transcripts_DBSCAN(cluster, 1);
		collapse_final_transcripts(cluster, 1);
	}

	// Assign Transcripts to Genes
	void assign_transcripts(std::shared_ptr<ClusterList> &cluster) {
		assign_to_genes(annotation, cluster, this -> get_contig_name(), 0); // Forward
		assign_to_genes(annotation, cluster, this -> get_contig_name(), 1); // Forward
	}

	/////////////////////////////////////////////////////////////
	/* Output Functions */

	// Print Clusters as GTF
	void write_gtf(std::ofstream &gtfFile) {
		if (ignore) { return; }
		list_ptr -> write_clusters_as_GTF(gtfFile);
		list_ptr.reset();
	}

	void get_stats(std::shared_ptr<ClusterList> &cluster) {
		assigned_reads = cluster -> get_assigned_reads();
		unassigned_reads = cluster -> get_unassigned_reads();
		ambiguous_reads = cluster -> get_ambiguous_reads();
		multimapped_reads = cluster -> get_multimapped_reads();
		low_quality_reads = cluster -> get_low_quality_reads();
		total_reads = cluster -> get_total_reads();
		transcript_num = cluster -> get_transcript_num();
	}

	/////////////////////////////////////////////////////////////
	/* Thread Launcher */

	void launch() {
		this -> open_alignment_file();                   // open files
		{
			std::shared_ptr<ClusterList> cluster;
			this -> create_clusters(cluster);            // find clusters
			this -> close_alignment_file();              // close files
			if (!ignore) {
				this -> collapse_clusters(cluster);      // collapse clusters
				this -> find_transcripts(cluster);       // dbscan clustering algorithm
				this -> assign_transcripts(cluster);     // overlap genes
				this -> list_ptr = cluster;              // Save for later
			}
			this -> get_stats(cluster);                  // Get Read Stats
		}
	}
};

