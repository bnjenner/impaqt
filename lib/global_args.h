///////////////////////////////////////
// Arguments Data Structure
#ifndef IMPAQT_ARGUMENTS_H
#define IMPAQT_ARGUMENTS_H

/*
    For the record, having this as a globally accessible struct MAY impact performance (albiet neglibible)
        for the following reasons.

    1.) Not hot on the cache, as it's infrequently accessed (slow relative to stack allocation)
    2.) ?


    I ultimately think splitting hairs over a few seconds is not necessarily optimal. If anything, readability
        and maintainability should be prioritized.
*/
namespace ImpaqtArguments {
struct GlobalArgs {

    // Files
    std::string alignment_file;         // sam or bam file
    std::string index_file;             // index filed
    std::string annotation_file;        // gff or gtf file

    // Program
    int threads;                        // threads

    // Alignments
    std::string library_type;           // library type (SE or PE)
    std::string stranded;               // strandedness
    bool nonunique_alignments;          // count primary and secondary alignments
    int mapq;                           // minimum mapq score
    int min_count;                      // min read count for dbscan
    int count_percentage;                // read count percentage for core reads
    int epsilon;                        // epsilon for dbscan

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