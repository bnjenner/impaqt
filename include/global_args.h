///////////////////////////////////////
// Arguments Data Structure
#ifndef IMPAQT_ARGUMENTS_H
#define IMPAQT_ARGUMENTS_H

namespace ImpaqtArguments {
struct GlobalArgs {

    // Files
    std::string alignment_file;         // sam or bam file
    std::string index_file;             // index filed
    std::string annotation_file;        // gff or gtf file

    // Program
    int threads;                        // threads

    // Alignments
    std::string stranded;               // strandedness
    bool nonunique_alignments;          // count primary and secondary alignments
    int mapq;                           // minimum mapq score
    int window_size;                    // window size
    int min_count;                      // min read count for dbscan
    int count_percentage;               // read count percentage for core reads
    int epsilon;                        // epsilon for dbscan
    double density_threshold;           // read density threshold for identifying organelle genomes

    // Features
    bool isGFF = false;                 // annotaiton file is gff
    std::string feature_tag;            // name of feature tag
    std::string feature_id;             // ID of feature

    // Output
    std::string gtf_output;             // name of output gtf file
};

extern GlobalArgs Args;
}
#endif
