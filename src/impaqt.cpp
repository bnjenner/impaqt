#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <vector>
#include <thread>
#include <condition_variable>
#include <mutex>
#include "api/BamReader.h"
#include "api/BamAux.h"
#include "parser.h"
#include "node.h"
#include "annotation.h"
#include "cluster.h"
#include "impaqt.h"
#include "queue.h"

// Mutex and Conditinoal vars
std::mutex main_mut;
std::condition_variable main_cv;
bool MAIN_THREAD = false;

//////////////////////////////////////
// Main
int main(int argc, char const ** argv) {

    auto start = std::chrono::high_resolution_clock::now();

    // Welcome!
    std::cerr << "[IMPAQT]\n";

    // Parse arguments
    ImpaqtArguments args;
    seqan::ArgumentParser::ParseResult res = argparse(argc, argv, args);

    if (res != seqan::ArgumentParser::PARSE_OK) {
        return res;
    }


    // Begin Parsing Files
    std::cerr << "[Parsing Input Files...]\n";

    std::cerr << "[...Annotation File...]\n";
    AnnotationList init_annotation(&args);
    // init_annotation.create_gene_graph();

    std::cerr << "[...Alignment File...]\n";
    Impaqt init_process(&args, 0);
    init_process.open_alignment_file();
    init_process.set_chrom_order();
    init_process.close_alignment_file();


    // Number of contigs for subdividing work across multiple threads
    // const int n = init_process.contig_map.size();
    const int n = 1; // only do first chromosome

    // Multithreading init
    std::vector<Impaqt*> processes;
    processes.reserve(n);

    for (int i = 0; i < n; i++) {
        processes.emplace_back(new Impaqt(&args, i));
        processes[i] -> copy_order(init_process.contig_map, init_process.contig_lengths);
        processes[i] -> copy_annotation(init_annotation, i);
    }


    // Send it
    std::cerr << "[Processing Data...]\n";

    int i = 0;
    const int proc = std::max(args.threads - 1, 1);
    {
        // initialize dispatch queue with N threads
        thread_queue call_queue(proc);
        do {
            // populate dispatch queue with necessary jobs
            while (i < n) {
                call_queue.dispatch([&, i] {processes[i] -> launch();}); // dispatch job
                i++;
            }

            // Wait for queue to be emptied
        } while (!call_queue.finished());
    }

    // std::unique_lock<std::mutex> main_lock(main_mut);   // lock main thread
    // main_cv.wait(main_lock, [] {return MAIN_THREAD;});  // wait for thread_queue destructor to let go
    // main_lock.unlock();                                 // unlock thread


    // // Summary Statistics
    // size_t total_ambiguous = 0;
    // size_t total_multimapping = 0;
    // size_t total_no_feature = 0;
    // size_t total_low_quality = 0;
    // size_t total_unique = 0;
    // size_t total_reads = 0;


    // // Write Results
    // std::cerr << "[Writing Results...]\n";
    // std::cerr << "[...Counts Data...]\n";


    // for (int i = 0; i < n; i++) {
    //     processes[i] -> print_counts();
    //     total_ambiguous += processes[i] -> ambiguous_reads;
    //     total_unique += processes[i] -> unique_reads;
    //     total_multimapping += processes[i] -> multimapped_reads;
    //     total_no_feature += processes[i] -> unassigned_reads;
    //     total_reads += processes[i] -> total_reads;
    // }

    // std::cout << "__unique\t" << total_unique << "\n";
    // std::cout << "__ambiguous\t" << total_ambiguous << "\n";
    // std::cout << "__multimapping\t" << total_multimapping << "\n";
    // std::cout << "__unassigned\t" << total_no_feature << "\n";
    // std::cout << "__total\t" << total_reads << std::endl;

    // // Output read cluster gtf if specified
    // if (args.gtf_output != "") {
    //     std::cerr << "[...Output GTFs...]\n";
    //     std::ofstream newFile;
    //     newFile.open(args.gtf_output);
    //     for (const auto &p : processes) { p -> print_gtf(); }
    //     newFile.close();
    // }


    // The longest line of "get the time" I have ever seen.
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

    // Say goodbye :)
    std::cerr << "[Program Complete!]\n";
    std::cerr << "[Runtime: " << duration.count() << " seconds]" << std::endl;

    return 0;
}
