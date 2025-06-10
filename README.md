# IMPAQT
# IN DEVELOPMENT
Currently only performs transcript identification and naive quantification. Gene assignment and more sophisticated quantification methods remain on the TO-DO list.

## Introduction

IMPAQT (Identifies Multiple Peaks and Quantifies Transcripts) is an
expression quantification method for TAGseq experiments developed 
by Bradley Jenner for his Undergraduate Honors Thesis at UC Davis. 
It operates on assumptions made about the distribution of sequencing reads
along the 3' UTR of a gene. Clustering these reads enables identification 
and quantification of gene level and transcript level expression for 
isoforms with distinct terminal exon usage. This method is particularly
useful in non-model organisms where 3' UTRs for most genes are poorly
annotated, resulting in massive data loss. It also can generate a GTF file
defining the boundaries and expression levels for each identified cluster. 
It relies partially on the bamtools and seqan C++ libraries, although mainly 
for argument parsing and reading in bam/gtf files. 

Additionally, the threadsafe queue and DBSCAN implementations were made possible 
due to the great work of github users EmbeddedArtistry and user Eleobert. 

For questions or comments, please contact
Bradley Jenner at <bnjenner@ucdavis.edu>

## Installation

0. Make sure cmake and make are installed on your machine.

1. Clone this repository and change into it.
```
git clone https://github.com/bnjenner/impaqt.git
cd impaqt
```

2. Create a build directory and change into it.
```
mkdir build
cd build
```

3. Compile
```
cmake ../
make
```

4. Add path to bash profile
```
echo "export PATH=$PATH:path/to/build_directory" >> ~/.bash_profile
source ~/.bash_profile
```
5. Give it a go! 

## Usage
```
SYNOPSIS
    impaqt input.sorted.bam annotation.gtf [options]

DESCRIPTION
    Identifies Multiple Peaks and Qauntifies Transcripts. Identifies and quantifies isoforms
    utilizing distinct terminal exons. Generates a counts file written to stdout and optionally a GTF
    file of identified read clusters.

REQUIRED ARGUMENTS
    BAM INPUT_FILE
    GTF INPUT_FILE

OPTIONS
    -h, --help
          Display the help message.
    -t, --threads INTEGER
          Number of processers for multithreading. Default: 1.
    -l, --library-type STRING
          Library type. Paired end is not recommended. Only used to check proper pairing. One of single and paired. Default: single.
    -s, --strandedness STRING
          Strandedness of library. One of forward and reverse. Default: forward.
    -n, --nonunique-alignments
          Count primary and secondary read alignments.
    -q, --mapq-min INTEGER
          Minimum mapping quality score to consider for counts. Default: 1.
    -m, --min-count INTEGER
          Minimum read count for DBSCAN transcript identification algorithm. (Minimum of 10) Default: 10.
    -p, --count-percentage INTEGER
          Minimum read count percentage for identifying core reads in DBSCAN algorithm. This will be the threshold unless percentage is less than 20. Default: 20.
    -e, --epsilon INTEGER
          Distance (in base pairs) for DBSCAN algorithm. Default: 150.
    -f, --feature-tag STRING
          Name of feature tag. Default: exon.
    -i, --feature-id STRING
          ID of feature (use for GFFs). Default: gene_id.
    -o, --output-gtf STRING
          Output read cluster GTF file and specify name.
    --version
          Display version information.

VERSION
    Last update: June 2025
    impaqt version: dev0
    SeqAn version: 2.4.0
```
