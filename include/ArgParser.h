#pragma once

#include <iostream>
#include <string>
#include <stdexcept>

#include "global_args.h"
#include "utils.h"

/////////////////////////////////////////////////////////////
// Argument Parser (hand-rolled; replaces the former seqan dependency)

// Result of parsing: proceed, exit with failure, or exit cleanly (help/version).
enum class ParseStatus { Ok, Error, Done };

namespace argparse_detail {

// Extension after the final '.' in the basename (no leading dot); "" if none.
// Matches the seqan getFileExtension behaviour the checks below relied on.
inline std::string file_extension(const std::string &path) {
    const size_t slash = path.find_last_of("/\\");
    const size_t dot = path.find_last_of('.');
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) { return ""; }
    return path.substr(dot + 1);
}

// Parse an integer option value; reports and returns false on bad input.
inline bool parse_int(const std::string &s, int &out, const std::string &name) {
    try {
        size_t pos = 0;
        const int v = std::stoi(s, &pos);
        if (pos != s.size()) { throw std::invalid_argument("trailing characters"); }
        out = v;
        return true;
    } catch (const std::exception &) {
        std::cerr << "ERROR: Option \"" << name << "\" expects an integer, got \"" << s << "\".\n";
        return false;
    }
}

// Parse a double option value; reports and returns false on bad input.
inline bool parse_double(const std::string &s, double &out, const std::string &name) {
    try {
        size_t pos = 0;
        const double v = std::stod(s, &pos);
        if (pos != s.size()) { throw std::invalid_argument("trailing characters"); }
        out = v;
        return true;
    } catch (const std::exception &) {
        std::cerr << "ERROR: Option \"" << name << "\" expects a number, got \"" << s << "\".\n";
        return false;
    }
}

inline void print_usage() {
    std::cerr <<
        "impaqt -- Identifies Multiple Peaks and Quantifies Transcripts.\n"
        "Identifies and quantifies isoforms utilizing distinct 3' ends. Generates a GTF\n"
        "file of identified transcripts and optionally a counts file to stdout if a\n"
        "reference annotation is provided.\n\n"
        "Usage: impaqt input.sorted.bam [options]\n\n"
        "Options:\n"
        "  -t, --threads INT             Number of processers for multithreading. [1]\n"
        "  -a, --annotation FILE         Annotation file (GTF or GFF). If set, a counts\n"
        "                                table is written to stdout. Type from extension. []\n"
        "  -s, --strandedness STR        Strandedness of library: forward or reverse. [forward]\n"
        "  -n, --nonunique-alignments    Count primary and secondary read alignments.\n"
        "  -q, --mapq-min INT            Minimum mapping quality score to consider. [1]\n"
        "  -w, --window-size INT         Window size to partition genome for read collection. [1000]\n"
        "  -m, --min-count INT           Min read count to initiate DBSCAN. (Hard minimum 10) [25]\n"
        "  -p, --count-percentage INT    Min read count percentage for core reads in DBSCAN. [5]\n"
        "  -e, --epsilon INT             Neighbor distance (bp) for DBSCAN. [50]\n"
        "  -d, --density-threshold DBL   Read density (#reads/#bp) to skip identification. [0]\n"
        "  -f, --feature-tag STR         Name of feature in GTF for assignment. [exon]\n"
        "  -u, --utr-tag STR             Name of UTR feature in GTF for assignment. [UTR]\n"
        "  -i, --feature-id STR          ID of feature to use for assignment. [gene_id]\n"
        "  -o, --output-gtf STR          Output GTF name. [BAM name + \".gtf\"]\n"
        "  -h, --help                    Print this help message and exit.\n"
        "      --version                 Print version and exit.\n";
}

} // namespace argparse_detail

inline ParseStatus argparse(int argc, char const **argv) {

    using namespace argparse_detail;

    // Defaults (formerly seqan setDefaultValue)
    ImpaqtArguments::Args.threads = 1;
    ImpaqtArguments::Args.annotation_file = "";
    ImpaqtArguments::Args.stranded = "forward";
    ImpaqtArguments::Args.nonunique_alignments = false;
    ImpaqtArguments::Args.mapq = 1;
    ImpaqtArguments::Args.window_size = 1000;
    ImpaqtArguments::Args.min_count = 25;
    ImpaqtArguments::Args.count_percentage = 5;
    ImpaqtArguments::Args.epsilon = 50;
    ImpaqtArguments::Args.density_threshold = 0;
    ImpaqtArguments::Args.feature_tag = "exon";
    ImpaqtArguments::Args.utr_tag = "UTR";
    ImpaqtArguments::Args.feature_id = "gene_id";
    ImpaqtArguments::Args.gtf_output = "";
    ImpaqtArguments::Args.isGFF = false;

    std::string bam;
    bool have_bam = false;

    for (int i = 1; i < argc; i++) {

        std::string tok = argv[i];

        if (tok == "-h" || tok == "--help") { print_usage(); return ParseStatus::Done; }
        if (tok == "--version") { std::cerr << "impaqt beta\n"; return ParseStatus::Done; }

        // Boolean flag (no value)
        if (tok == "-n" || tok == "--nonunique-alignments") {
            ImpaqtArguments::Args.nonunique_alignments = true;
            continue;
        }

        // Positional argument (the input BAM)
        if (tok.empty() || tok[0] != '-') {
            if (have_bam) {
                std::cerr << "ERROR: Unexpected extra argument: \"" << tok << "\".\n";
                return ParseStatus::Error;
            }
            bam = tok;
            have_bam = true;
            continue;
        }

        // Option with a value. Support both "--name value" and "--name=value".
        std::string name = tok, inline_val;
        bool has_inline = false;
        const size_t eq = tok.find('=');
        if (tok.rfind("--", 0) == 0 && eq != std::string::npos) {
            name = tok.substr(0, eq);
            inline_val = tok.substr(eq + 1);
            has_inline = true;
        }

        // Pull the value from "=" or the next token.
        auto get_value = [&](std::string &out) -> bool {
            if (has_inline) { out = inline_val; return true; }
            if (i + 1 >= argc) {
                std::cerr << "ERROR: Option \"" << name << "\" requires a value.\n";
                return false;
            }
            out = argv[++i];
            return true;
        };

        std::string val;
        if (name == "-t" || name == "--threads") {
            if (!get_value(val) || !parse_int(val, ImpaqtArguments::Args.threads, name)) { return ParseStatus::Error; }

        } else if (name == "-a" || name == "--annotation") {
            if (!get_value(ImpaqtArguments::Args.annotation_file)) { return ParseStatus::Error; }

        } else if (name == "-s" || name == "--strandedness") {
            if (!get_value(ImpaqtArguments::Args.stranded)) { return ParseStatus::Error; }
            if (ImpaqtArguments::Args.stranded != "forward" && ImpaqtArguments::Args.stranded != "reverse") {
                std::cerr << "ERROR: --strandedness must be \"forward\" or \"reverse\".\n";
                return ParseStatus::Error;
            }

        } else if (name == "-q" || name == "--mapq-min") {
            if (!get_value(val) || !parse_int(val, ImpaqtArguments::Args.mapq, name)) { return ParseStatus::Error; }

        } else if (name == "-w" || name == "--window-size") {
            if (!get_value(val) || !parse_int(val, ImpaqtArguments::Args.window_size, name)) { return ParseStatus::Error; }

        } else if (name == "-m" || name == "--min-count") {
            if (!get_value(val) || !parse_int(val, ImpaqtArguments::Args.min_count, name)) { return ParseStatus::Error; }

        } else if (name == "-p" || name == "--count-percentage") {
            if (!get_value(val) || !parse_int(val, ImpaqtArguments::Args.count_percentage, name)) { return ParseStatus::Error; }

        } else if (name == "-e" || name == "--epsilon") {
            if (!get_value(val) || !parse_int(val, ImpaqtArguments::Args.epsilon, name)) { return ParseStatus::Error; }

        } else if (name == "-d" || name == "--density-threshold") {
            if (!get_value(val) || !parse_double(val, ImpaqtArguments::Args.density_threshold, name)) { return ParseStatus::Error; }

        } else if (name == "-f" || name == "--feature-tag") {
            if (!get_value(ImpaqtArguments::Args.feature_tag)) { return ParseStatus::Error; }

        } else if (name == "-u" || name == "--utr-tag") {
            if (!get_value(ImpaqtArguments::Args.utr_tag)) { return ParseStatus::Error; }

        } else if (name == "-i" || name == "--feature-id") {
            if (!get_value(ImpaqtArguments::Args.feature_id)) { return ParseStatus::Error; }

        } else if (name == "-o" || name == "--output-gtf") {
            if (!get_value(ImpaqtArguments::Args.gtf_output)) { return ParseStatus::Error; }

        } else {
            std::cerr << "ERROR: Unknown option \"" << name << "\".\n";
            return ParseStatus::Error;
        }
    }

    // Require the positional BAM
    if (!have_bam) {
        std::cerr << "ERROR: Missing required input BAM file.\n";
        print_usage();
        return ParseStatus::Error;
    }

    // Check file type of the input alignment
    if (file_extension(bam) != "bam") {
        std::cerr << "ERROR: Unaccepted File Format: \"." << file_extension(bam)
                  << "\". Only accepts \".bam\" extension.\n";
        return ParseStatus::Error;
    }

    ImpaqtArguments::Args.alignment_file = bam;
    ImpaqtArguments::Args.index_file = bam + ".bai";

    if (!file_exists(ImpaqtArguments::Args.alignment_file)) {
        std::cerr << "ERROR: Alignment file \"" << ImpaqtArguments::Args.alignment_file << "\" does not exist.\n";
        throw std::runtime_error("ERROR: Make sure alignment file exists.");
    }

    if (ImpaqtArguments::Args.gtf_output == "") {
        ImpaqtArguments::Args.gtf_output = ImpaqtArguments::Args.alignment_file + ".gtf";
    }

    // Check file type of the annotation
    if (ImpaqtArguments::Args.annotation_file != "") {
        if (!file_exists(ImpaqtArguments::Args.annotation_file)) {
            std::cerr << "ERROR: Annotation file \"" << ImpaqtArguments::Args.annotation_file << "\" does not exist.\n";
            throw std::runtime_error("ERROR: Make sure annotation file exists.");
        }
        const std::string ext = file_extension(ImpaqtArguments::Args.annotation_file);
        if (ext != "gtf" && ext != "gff") {
            std::cerr << "ERROR: Unaccepted File Format: \"." << ext
                      << "\". Only accepts \".gtf\" and \".gff\" extension.\n";
            return ParseStatus::Error;
        }
        if (ext == "gff") { ImpaqtArguments::Args.isGFF = true; }
    }

    return ParseStatus::Ok;
}
