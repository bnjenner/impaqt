#include "AnnotationList.h"
#include "ClusterList.h"
#include "utils.h"
#include "DBSCAN.h"
#include "AssignClusters.h"

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/* Impaqt Process Class */

class Impaqt {

private:

	// Alignment file Readers
	BamTools::BamReader inFile;						// Bam File Object
	BamTools::BamAlignment alignment;					// BamAlignmentRecord record;

	static AnnotationList annotation;					// Annotation list for genes
	static std::string alignment_file_name;					// alignment file
	static std::string index;						// alignment index file
	int chrom_index;							// chromosome number
	bool ignore_chr = false;						// to ignore for downstream

	static std::unordered_map<int, std::string> contig_map;			// Links Index to Contig Name
	static std::unordered_map<int, int> contig_lengths;			// Links Index to Contig Length

	ClusterList* cluster_list;						// List for clusters

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
	Impaqt(const int &chrom_index) {
		alignment_file_name = ImpaqtArguments::Args.alignment_file;
		index = ImpaqtArguments::Args.index_file;
		this -> chrom_index = chrom_index;
	}

	// Destructor
	~Impaqt() { if (!ignore_chr) { delete cluster_list; } }


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
	ClusterList* get_clusters() { return cluster_list; }

	// Get Chromosome Info
	int get_chrom_index() { return chrom_index; }
	int get_chrom_num() { return contig_map.size(); }
	bool is_ignored() { return ignore_chr; }
	std::string get_contig_name() { return contig_map[chrom_index]; }
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
		if (!inFile.OpenIndex(index)) {
			std::cerr << "ERROR: Could not read index file: " << index << "\n";
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
	void create_clusters() {

		cluster_list = new ClusterList();
		cluster_list -> initialize(chrom_index, contig_map[chrom_index], contig_lengths[chrom_index]);

		if (!inFile.Jump(chrom_index)) {
			std::cerr << "//ERROR: Could not jump to region: " << this -> get_contig_name() << "\n";
			throw "ERROR: Could not jump to region. Make sure BAM header is correct.";
		}

		// If failed to create clusters, flag to ignore
		if (!(cluster_list -> create_clusters(inFile, alignment))) { ignore_chr = true; }
	}

	// Merge neighboring clusters and remove zeroes
	void collapse_clusters() {
		cluster_list -> collapse_clusters(0); // Forward
		cluster_list -> collapse_clusters(1); // Reverse
	}

	// Differentiate Transcripts
	void find_transcripts() {
		find_transcripts_DBSCAN(cluster_list, 0); // Forward
		find_transcripts_DBSCAN(cluster_list, 1); // Reverse
	}

	// Assign Transcripts to Genes
	void assign_transcripts() {
		assign_to_genes(annotation, cluster_list, contig_map[chrom_index], 0); // Forward
		assign_to_genes(annotation, cluster_list, contig_map[chrom_index], 1); // Forward
	}


	/////////////////////////////////////////////////////////////
	/* Output Functions */

	// Print Clusters as GTF
	void write_gtf(std::ofstream &gtfFile) {
		if (ignore_chr) { return; }
		cluster_list -> write_clusters_as_GTF(gtfFile);
	}

	void get_stats() {
		assigned_reads = cluster_list -> get_assigned_reads();
		unassigned_reads = cluster_list -> get_unassigned_reads();
		ambiguous_reads = cluster_list -> get_ambiguous_reads();
		multimapped_reads = cluster_list -> get_multimapped_reads();
		low_quality_reads = cluster_list -> get_low_quality_reads();
		total_reads = cluster_list -> get_total_reads();
		transcript_num = cluster_list -> get_transcript_num();
	}


	/////////////////////////////////////////////////////////////
	/* Thread Launcher */

	void launch() {
		this -> open_alignment_file();				// open files
		this -> create_clusters();				// find clusters
		this -> close_alignment_file();			// close files
		if (!ignore_chr) {
			this -> collapse_clusters();			// collapse clusters
			this -> find_transcripts();			// dbscan clustering algorithm
			this -> assign_transcripts();			// overlap genes
		}
		this -> get_stats();					// Get Read Stats
		if (ignore_chr) { delete cluster_list; }
	}
};

