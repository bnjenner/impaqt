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
	BamTools::BamReader inFile;
	BamTools::BamAlignment alignment;

	// Files and Data Structure
	ClusterList* cluster_list;
	static AnnotationList annotation;
	static std::string alignment_file_name;
	static std::string index_file_name;
	static std::unordered_map<int, std::string> contig_map;
	static std::unordered_map<int, int> contig_lengths;
	
	std::string contig_name;
	int contig_index;
	int contig_length;
	bool ignore = false;

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

	// Destructor
	~Impaqt() {}

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
		if (!inFile.Open(alignment_file_name)) {
			std::cerr << "ERROR: Could not read alignment file: " << alignment_file_name << "\n";
			throw "ERROR: Make sure alignment file exists.";
		}
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

		// Generate Contig Map (contig indicies -> contig name)
		BamTools::RefVector references = inFile.GetReferenceData();
		for (int i = 0; i < references.size(); i++) {
			contig_map[i] = references.at(i).RefName;
			contig_lengths[i] = references.at(i).RefLength;
		}
	}

	void set_contigs() {
		contig_name = contig_map[contig_index];
		contig_length = contig_lengths[contig_index];
	}

	void add_annotation() {
		annotation = AnnotationList();
		annotation.create_gene_list();
	}

	/////////////////////////////////////////////////////////////
	/* Cluster Related Functions */

	void create_clusters() {

		cluster_list = new ClusterList(contig_index, contig_name, contig_length);

		if (!inFile.Jump(contig_index)) {
			std::cerr << "//ERROR: Could not jump to region: " << contig_name << "\n";
			throw "ERROR: Could not jump to region. Make sure BAM header is correct.";
		}

		// If failed to create clusters, flag to ignore
		if (!(cluster_list -> create_clusters(inFile, alignment))) { ignore = true; }
	}

	void collapse_clusters() {
		int t_strand = 0; // Forward
		cluster_list -> collapse_clusters(t_strand);
		cluster_list -> collapse_clusters(!t_strand);
	}

	void find_transcripts() {
		if (ignore) { return; }
		int t_strand = 0; // Forward
		identify_transcripts(cluster_list, t_strand);
		identify_transcripts(cluster_list, !t_strand);
	}

	void assign_transcripts() {
		int t_strand = 0; // Forward
		assign_to_genes(annotation, cluster_list, contig_name, t_strand);
		assign_to_genes(annotation, cluster_list, contig_name, !t_strand);
	}

	/////////////////////////////////////////////////////////////
	/* Output Functions */

	void write_gtf(std::ofstream &gtfFile) {
		if (ignore) { return; }
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
		this -> set_contigs();
		this -> open_alignment_file();
		this -> create_clusters();
		this -> close_alignment_file();
		if (!ignore) {
			this -> collapse_clusters();
			this -> find_transcripts();
			if (ImpaqtArguments::Args.annotation_file != "") {
				this -> assign_transcripts();
			}
		}
		this -> get_stats();
	}
};

