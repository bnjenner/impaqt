# Impaqt [![C/C++ CI](https://github.com/bnjenner/impaqt/actions/workflows/c-cpp.yml/badge.svg)](https://github.com/bnjenner/impaqt/actions/workflows/c-cpp.yml)

## Areas to Improve

### Parallelization
Current implemention of multithreading splits the workload by contigs. While this can help chew through
highly fragmented genomes, it does not improve performance on individual contigs that may have high
coverage. DBSCAN seems to be the bottlenck for this situtaion, as the O(n2) complexity of this algorithm
really chugs with these regions. Parallelizing seems doable. 

### Default Parameter Values
There are three main problem parameters currently: window-size, count-percentage, and epsilon. The trade off with
window-size is fairly straight forward, the smaller window size, the higher the resolution when identifying distinct
transcripts; however, a longer window increases false positives, meaning transcripts are considered separate when
in actuality they might not be. The inverse problem is true for the epsilon parameter, which defines the distance
threshold between points to be considered "neighbors" for core point identification.

The last parameter is trickier, because both identification and quantification of transcripts is dependent on the 
--count-percentage parameter. This is because this percentage is used to define the concentration parameter in DBSCAN 
that determines core points. Setting this number statically for a region that is dominated by a single transcript will
prevent identification of more lowly expressed neighboring transcripts, even when it makes sence biologically that this 
transcript is in fact unique (i.e. a transcript located downstream in the UTR). Quantification of transcripts is also impacted
because it occurs after and depends on the number of core points identified. Briefly, an isoform's counts value is equal to 
the proportion of the total core points contained in that transcript multiplied by the total read counts for that region. Luckily, 
this does not necessarily scale linearly with the length of the transcript since core points are identified using density, so
a pseudo-length bias is not an issue as with standard RNAseq. However, the core point concentration parameter could be biasing
identification of transcripts (and subsequently biasing quantification as that is dependent on the conc. parameter), so we
may want to explore more sophisticated ways of setting that parameter or instead settting it liberally only to refine the results
probabilistically after. 

### Transcript Bounds
Currently, transcript identification finds the bounds by considering all points assigned to a particular group by
DBSCAN. Usually, this is fine; however, noise in the form of spuriously aligned reads or unannotated exons result
in bounds that do not always align exonic starts/stops. This can somewhat hurt assignment, but usually in the case
of small RNAs located in intronic regions.

### Small RNA Species
For some reason, small RNA's love to show up in these datasets and cause problems, particularly in intronic regions. This algorithm
was designed to function independent of a reference annotation while also seeking to identify mRNA transcripts, not the extremely
narrow peaks resulting from small RNA expression. Impaqt will not see these reads as intronic or belonging to a small RNA, it will
treat them as if they were another exon, even though they do not belong to the same gene. This can sometimes result in erroneous
ambiguous assignment.

### Dependencies
Dependency on seqan is sort of not necessary. Bamtools is a necessary dependency, but I need to find a better 
way to include it besides just having an 'ext' folder containing the files I need. Additionally, the Bamtools
version is a bit behind and calls to the API will need to be revised. 

### Test Code Coverage
Need I say more.
