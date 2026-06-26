# Impaqt [![C/C++ CI](https://github.com/bnjenner/impaqt/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/bnjenner/impaqt/actions/workflows/c-cpp.yml)
Currently performs transcript identification, gene assignment, and naive quantification. 
A more sophisticated quantification method remain on the TO-DO list.

![ ](./docs/example.png)

## Introduction

IMPAQT (Identifies Multiple Peaks and Quantifies Transcripts) is a transcript
identification and gene expression quantification method for TAGseq and
3' mRNAseq experiments. It operates on assumptions about the distribution 
of reads along the 3' UTR of expressed genes. Clustering these reads 
enables pseudo-transcript identification and quantification of expression at the 
gene and transcript level for isoforms utilizing distinct 3' ends. 

It generates a GTF file defining the boundaries of each transcript and their 
expression level as well as, optionally, a gene expression counts table 
if a reference annotation is provided. 

This method is particularly useful in non-model organisms where 3' UTRs for 
most genes are poorly annotated (resulting in massive data loss). Increased
gene density also tends to hurt the assignment of transcripts by this 
aglorithm, as it increases assignment ambiguity. Reads for which a reasonable
transcript of origin cannot be identified are handled individually. 

## Installation

0. Make sure `cmake` (≥ 3.18), `make`, a C++17 compiler, and `zlib` are installed on your machine.

```
# Linux
sudo apt install cmake zlib1g-dev

# Mac
brew install cmake zlib
```

> **Note:** The remaining third-party libraries (bamtools, GoogleTest) are
> downloaded and built automatically by CMake on the first configure, so an
> internet connection is required the first time you build. They are pinned to
> specific versions and are not vendored in this repository.

1. Clone this repository and change into it.
```
git clone https://github.com/bnjenner/impaqt.git
cd impaqt
```

### Quick install (recommended)

Use the provided script, which configures, builds, and (optionally) installs:
```
./install.sh --install            # build + install to the system prefix (uses sudo)
./install.sh                      # build only -> ./build/impaqt
./install.sh --install --prefix ~/.local   # install without sudo
```
Run `./install.sh --help` for all options (build type, parallel jobs, prefix).

### Manual build

```
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
sudo cmake --install build         # installs only the impaqt binary
```

Then give it a go!
```
impaqt input.sorted.bam
```

## Usage
```
SYNOPSIS
    impaqt input.sorted.bam [options]

DESCRIPTION
    Identifies Multiple Peaks and Qauntifies Transcripts. Identifies and quantifies isoforms utilizing distinct 3'
    ends. Generates a GTF file of identified transcripts and optionally a counts file written to stdout if a reference
    annotation is provided.

REQUIRED ARGUMENTS
    BAM INPUT_FILE

OPTIONS
    -h, --help
          Display the help message.
    -t, --threads INTEGER
          Number of processers for multithreading. Default: 1.
    -a, --annotation INPUT_FILE
          Annotation File (GTF or GFF). If specified, a counts table will be output through standard out. NOTICE: File
          type identified by file extension. Default: .
    -s, --strandedness STRING
          Strandedness of library. One of forward and reverse. Default: forward.
    -n, --nonunique-alignments
          Count primary and secondary read alignments.
    -q, --mapq-min INTEGER
          Minimum mapping quality score to consider. Default: 1.
    -w, --window-size INTEGER
          Window size to use to parition genome for read collection. Default: 1000.
    -m, --min-count INTEGER
          Minimum read count to initiate DBSCAN transcript identification algorithm. (Hard minimum of 10) Default: 25.
    -p, --count-percentage INTEGER
          Minimum read count percentage for identifying core reads in DBSCAN algorithm. This will be the threshold
          unless number of reads is less than 10. Default: 5.
    -e, --epsilon INTEGER
          Distance (in base pairs) for neighboring reads in DBSCAN algorithm. This should generally be 0.5-1.5x the
          read length, depending on desired isoform sensitivity (lower = more sensitive). Default: 50.
    -d, --density-threshold DOUBLE
          Read density threshold (# reads / # bps) to skip transcript identification. Assignment in super dense
          regions (usually the mitochrondria) doesn't really benefit from transcript identificaiton. Default is unset.
          Default: 0.
    -f, --feature-tag STRING
          Name of feature in GTF for assignment. Default: exon.
     -u, --utr-tag STRING
          Name of UTR feature in GTF for assignment. Default: UTR.
    -i, --feature-id STRING
          ID of feature to use for feature assignment. Default: gene_id.
    -o, --output-gtf STRING
          Specify name of cluster GTF file. Default is BAM name + ".gtf".
    --version
          Display version information.

VERSION
    Last update: August 2025
    impaqt version: beta
```

## Dependencies
Build requirements: `cmake` (≥ 3.18), a C++17 compiler, `make`, and `zlib`.

The following libraries are fetched and built automatically by CMake (via
`FetchContent`) on the first configure — they are pinned to specific versions and
not vendored in this repository:
- [bamtools](https://github.com/pezmaster31/bamtools) (v2.5.3) — BAM file I/O.
- [GoogleTest](https://github.com/google/googletest) (release-1.12.1) — unit tests only.

The DBSCAN clustering algorithm is inspired by github user [Eleobert](https://github.com/Eleobert/dbscan/blob/master/dbscan.cpp), and the thread dispatch by [EmbeddedArtistry](https://github.com/embeddedartistry/embedded-resources/blob/master/examples/cpp/dispatch.cpp).

## Contact
For questions or comments, please contact
Bradley Jenner at <bnjenner@bu.edu>


