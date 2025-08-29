#include "seqan/arg_parse.h"

/////////////////////////////////////////////////////////////
// Argument Parser
seqan::ArgumentParser::ParseResult argparse(int argc, char const **argv) {

    // Setup ArgumentParser.
    seqan::ArgumentParser parser("impaqt");
    seqan::addDescription(parser,
                          "Identifies Multiple Peaks and Qauntifies Transcripts. Identifies and quantifies isoforms utilizing distinct 3' ends. Generates a GTF file of identified transcripts and optionally a counts file written to stdout if a reference annotation is provided.");


    // Define Arguments
    seqan::addArgument(parser, seqan::ArgParseArgument(
                           seqan::ArgParseArgument::INPUT_FILE, "BAM"));

    // Define Program Options
    seqan::addOption(parser, seqan::ArgParseOption(
                         "t", "threads",
                         "Number of processers for multithreading.",
                         seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::setDefaultValue(parser, "threads", "1");

    seqan::addOption(parser, seqan::ArgParseOption(
                         "a", "annotation",
                         "Annotation File (GTF or GFF). If specified, a counts table will be output through standard out. NOTICE: File type identified by file extension.",
                         seqan::ArgParseArgument::INPUT_FILE, "STRING"));
    seqan::setDefaultValue(parser, "annotation", "");


    // Define Read Options
    seqan::addOption(parser, seqan::ArgParseOption(
                         "s", "strandedness", "Strandedness of library.",
                         seqan::ArgParseArgument::STRING, "STRING"));
    seqan::setDefaultValue(parser, "strandedness", "forward");
    seqan::setValidValues(parser, "strandedness", "forward reverse");

    seqan::addOption(parser, seqan::ArgParseOption(
                         "n", "nonunique-alignments", "Count primary and secondary read alignments."));

    seqan::addOption(parser, seqan::ArgParseOption(
                         "q", "mapq-min",
                         "Minimum mapping quality score to consider.",
                         seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::setDefaultValue(parser, "mapq-min", "1");

    // Binning Options
    seqan::addOption(parser, seqan::ArgParseOption(
                  "w", "window-size",
                  "Window size to use to parition genome for read collection.",
                  seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::setDefaultValue(parser, "window-size", "1000");

    // Define DBSCAN Options
    seqan::addOption(parser, seqan::ArgParseOption(
                  "m", "min-count",
                  "Minimum read count to initiate DBSCAN transcript identification algorithm. (Hard minimum of 10)",
                  seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::setDefaultValue(parser, "min-count", "25");

    seqan::addOption(parser, seqan::ArgParseOption(
                  "p", "count-percentage",
                  "Minimum read count percentage for identifying core reads in DBSCAN algorithm. This will be the threshold unless number of reads is less than 10.",
                  seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::setDefaultValue(parser, "count-percentage", "5");

    seqan::addOption(parser, seqan::ArgParseOption(
                  "e", "epsilon",
                  "Distance (in base pairs) for neighboring reads in DBSCAN algorithm. This should generally be 0.5-1.5x the read length, depending on desired isoform sensitivity (lower = more sensitive).",
                  seqan::ArgParseArgument::INTEGER, "INT"));
    seqan::setDefaultValue(parser, "epsilon", "50");

    seqan::addOption(parser, seqan::ArgParseOption(
                  "d", "density-threshold",
                  "Read density threshold (# reads / # bps) to skip transcript identification. Assignment in super dense regions (usually the mitochrondria) doesn't really benefit from transcript identificaiton. Default is unset.",
                  seqan::ArgParseArgument::DOUBLE, "DOUBLE"));
    seqan::setDefaultValue(parser, "density-threshold", "0");


    // Define Annotation Options
    seqan::addOption(parser, seqan::ArgParseOption(
                         "f", "feature-tag", "Name of feature in GTF for assignment.",
                         seqan::ArgParseArgument::STRING, "STRING"));
    seqan::setDefaultValue(parser, "feature-tag", "exon");

    seqan::addOption(parser, seqan::ArgParseOption(
                         "i", "feature-id", "ID of feature to use for feature assignment.",
                         seqan::ArgParseArgument::STRING, "STRING"));
    seqan::setDefaultValue(parser, "feature-id", "gene_id");

    seqan::addOption(parser, seqan::ArgParseOption(
                         "o", "output-gtf", "Specify name of cluster GTF file. Default is BAM name + \".gtf\".",
                         seqan::ArgParseArgument::STRING, "STRING"));

    seqan::addUsageLine(parser, "input.sorted.bam [options]");
    seqan::setDefaultValue(parser, "version-check", "OFF");
    seqan::hideOption(parser, "version-check");
    seqan::setVersion(parser, "beta");
    seqan::setDate(parser, "August 2025");

    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    // Check if Parse was successful
    if (res != seqan::ArgumentParser::PARSE_OK) { return seqan::ArgumentParser::PARSE_ERROR; }

    // Check file type of first positional arg
    std::string input_file_ext = seqan::getFileExtension(getArgument(parser, 0));
    if (input_file_ext != "bam") {
        std::cerr << "ERROR: Unaccapetd File Format: \"." << input_file_ext <<  "\". Only accepts \".bam\",  extension.\n";
        return seqan::ArgumentParser::PARSE_ERROR;
    }

    // Get arguments
    seqan::getArgumentValue(ImpaqtArguments::Args.alignment_file, parser, 0);
    seqan::getArgumentValue(ImpaqtArguments::Args.index_file, parser, 0);
    ImpaqtArguments::Args.index_file = ImpaqtArguments::Args.index_file + ".bai";

    // Populate options
    seqan::getOptionValue(ImpaqtArguments::Args.annotation_file, parser, "annotation");
    seqan::getOptionValue(ImpaqtArguments::Args.threads, parser, "threads");
    seqan::getOptionValue(ImpaqtArguments::Args.stranded, parser, "strandedness");
    ImpaqtArguments::Args.nonunique_alignments = seqan::isSet(parser, "nonunique-alignments");
    seqan::getOptionValue(ImpaqtArguments::Args.mapq, parser, "mapq-min");
    seqan::getOptionValue(ImpaqtArguments::Args.window_size, parser, "window-size");
    seqan::getOptionValue(ImpaqtArguments::Args.min_count, parser, "min-count");
    seqan::getOptionValue(ImpaqtArguments::Args.count_percentage, parser, "count-percentage");
    seqan::getOptionValue(ImpaqtArguments::Args.epsilon, parser, "epsilon");
    seqan::getOptionValue(ImpaqtArguments::Args.density_threshold, parser, "density-threshold");
    seqan::getOptionValue(ImpaqtArguments::Args.feature_tag, parser, "feature-tag");
    seqan::getOptionValue(ImpaqtArguments::Args.feature_id, parser, "feature-id");
    seqan::getOptionValue(ImpaqtArguments::Args.gtf_output, parser, "output-gtf");


    if (ImpaqtArguments::Args.gtf_output == "") {
      ImpaqtArguments::Args.gtf_output = ImpaqtArguments::Args.alignment_file + ".gtf";
    }

    // Check file type of annotation
    if (ImpaqtArguments::Args.annotation_file != "") {
      input_file_ext = seqan::getFileExtension(getOption(parser, "annotation"));
      if (input_file_ext != "gtf" && input_file_ext != "gff") {
          std::cerr << "ERROR: Unaccapetd File Format: \"." << input_file_ext <<  "\". Only accepts \".gtf\" and \".gff\",  extension.\n";
          return seqan::ArgumentParser::PARSE_ERROR;
      }
      if (input_file_ext == "gff") { ImpaqtArguments::Args.isGFF = true; }
    }

    return seqan::ArgumentParser::PARSE_OK;
}
