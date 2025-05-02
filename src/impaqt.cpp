#include <iostream>
#include <fstream>
#include <chrono>
#include <string>
#include <vector>
#include <thread>
#include <condition_variable>
#include <mutex>
#include <global_args.h>
#include "parser.h"
#include "annotation.h"
#include "impaqt.h"
#include "queue.h"


// Globals
ImpaqtArguments::GlobalArgs ImpaqtArguments::Args;

// Static Member Defintions
std::string Impaqt::alignment_file_name;
std::string Impaqt::index;
std::unordered_map<int, std::string> Impaqt::contig_map;
std::unordered_map<int, int> Impaqt::contig_lengths;

// Mutex and Conditinoal vars
std::mutex main_mut;
std::condition_variable main_cv;
bool MAIN_THREAD = false;


//////////////////////////////////////
// Main
int main(int argc, char const ** argv) {

    auto start = std::chrono::high_resolution_clock::now();

    // Parse arguments
    seqan::ArgumentParser::ParseResult res = argparse(argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK) { return res; }

    // Welcome!
    std::cerr << "// IMPAQT\n";
    std::cerr << "// Parsing Input Files...\n";


    std::cerr << "//     Annotation File...\n";
    AnnotationList annotation; // Initiliazed using global args object 
    annotation.create_gene_graph();


    std::cerr << "//     Alignment File....\n";
    std::vector<Impaqt*> processes;
    processes.emplace_back(new Impaqt(0));
    processes[0] -> open_alignment_file();
    processes[0] -> set_chrom_order();

    // Number of contigs for subdividing work across multiple threads
    // int n = 1;
    int n = processes[0] -> get_chrom_num();

    // Multithreading init
    if (n > 1) {
        processes.reserve(n);
        for (int i = 1; i < n; i++) {
            processes.emplace_back(new Impaqt(i));
            // processes[i] -> copy_annotation(init_annotation, i);
        }
    }



    std::cerr << "// Processing Data.......\n";
    int i = 0;
    const int proc = std::max(ImpaqtArguments::Args.threads - 1, 1);
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

    std::unique_lock<std::mutex> main_lock(main_mut);   // lock main thread
    main_cv.wait(main_lock, [] {return MAIN_THREAD;});  // wait for thread_queue destructor to let go
    main_lock.unlock();                                 // unlock thread


    // // Summary Statistics
    // size_t total_ambiguous = 0;
    size_t total_multimapping = 0;
    // size_t total_no_feature = 0;
    // size_t total_low_quality = 0;
    // size_t total_unique = 0;
    size_t total_reads = 0;


    // Write Results
    std::cerr << "// Writing Results.......\n";
    std::cerr << "//     Counts Data.......\n";


    for (int i = 0; i < n; i++) {
        // processes[i] -> print_counts();
        // total_ambiguous += processes[i] -> ambiguous_reads;
        // total_unique += processes[i] -> unique_reads;
        total_multimapping += processes[i] -> get_multimapped_reads();
        // total_no_feature += processes[i] -> unassigned_reads;
        total_reads += processes[i] -> get_total_reads();
        delete processes[i];
    }

    // std::cout << "__unique\t" << total_unique << "\n";
    // std::cout << "__ambiguous\t" << total_ambiguous << "\n";
    std::cout << "// multimapping\t" << total_multimapping << "\n";
    // std::cout << "__unassigned\t" << total_no_feature << "\n";
    std::cout << "// total\t" << total_reads << std::endl;

    // // Output read cluster gtf if specifiedq
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
    std::cerr << "// Program Complete!\n";
    std::cerr << "// Runtime: " << duration.count() << " seconds" << std::endl;

    return 0;
}
