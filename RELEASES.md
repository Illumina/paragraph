# Paragraph Release Notes / Change Log

# Version 2.1

| Date Y-m-d | Ticket  | Description                                                          |
|------------|---------|----------------------------------------------------------------------|
| 2018-12-06 | GT-675  | Fix filters and alignment stats. Change depth test threshold on lower end |
| 2018-11-08 | GT-660  | Optimize GQ for variant genotypes                                    |
| 2018-11-02 | GT-656  | Improvement for simple SV genotyping                                 |
| 2018-07-19 | GT-501  | Breakpoint depth test based on normal distribution                   |
| 2018-07-16 | GT-539  | VCF now output genotypes for all samples in manifest and input VCF   |
| 2018-06-28 | GT-527  | --graph-sequence-matching yes fails with boost 1.63                  |

# Version 2.0

| Date Y-m-d | Ticket  | Description                                                          |
|------------|---------|----------------------------------------------------------------------|
| 2018-06-27 | GT-490  | Paragraph 2.0 release; disable Poisson depth test by default         |
| 2018-06-27 | GT-495  | Improved output of phasing information and paths                     |
| 2018-06-26 | GT-402  | support genotyping on male chrX                                      |
| 2018-06-07 | GT-496  | -A still causes grmpy to run out of memory                           |
| 2018-06-06 | GT-492  | -A causes grmpy to run out of memory                                 |
| 2018-05-31 | GT-491  | paragraph fails to compile with boost 1.53                           |
| 2018-05-31 | GT-484  | Output VCF with new genotypes as computed by grmpy                   |
| 2018-05-31 | GT-477  | Add features to output alignment in multigrmpy and perform superlocus VCF splitting |
| 2018-05-30 | GT-485  | Running GiaB VCF with klib alignment produces error message          |
| 2018-05-24 | GT-473  | Round floating point numbers in json output                          |
| 2018-05-24 | GT-471  | Report read-counts for path families instead of paths                |
| 2018-05-24 | GT-476  | some command line parameters ignored in paragraph                    |
| 2018-05-23 | GT-474  | support for command line location of graphtools source tree          |
| 2018-05-16 | GT-436  | post-process bad KmerAlignment results with smith waterman           |
| 2018-05-16 | GT-451  | Add path family based phasing and haplotype reconstruction           |
| 2018-05-15 | GT-445  | Add path aligner to produce faster alignments for exact matches      |
| 2018-05-14 | GT-462  | Minor fixes to vcf2paragraph.py for paralog graphs                   |
| 2018-05-10 | GT-463  | vcf2paragraph fails to produce alt paths                             |
| 2018-05-09 | GT-459  | Update to graph-tools 0.2.x                                          |
| 2018-04-23 | GT-450  | Local align Klib introducess odd soft clip                           |
| 2018-04-25 | GT-454  | harmonize output options between paragraph and grmpy                 |
| 2018-04-20 | GT-446  | Fix variant extraction due to bug in cigarToRefVar                   |
| 2018-04-20 | GT-449  | performance lost due to recent merge making bam/cram open for each graph instead of once per thead-file |
| 2018-04-19 | GT-447  | Fix multigrampy event ID                                             |
| 2018-04-16 | GT-439  | Fix logging + verbosity for multigrmpy.py                            |
| 2018-04-16 | GT-438  | Fix high memory usage in grmpy and paragraph                         |
| 2018-04-13 | GT-435  | Add support to vcf2paragraph.py for symbolic INS allele in SV VCFs   |
| 2018-04-13 | GT-437  | Fix multithreading race conditions                                   |
| 2018-04-11 | GT-353  | Refactor to use graph-tools library rather than protobuf             |
| 2018-04-07 | GT-417  | Refactor VCF to graph conversion                                     |

# Version 1.2

| Date Y-m-d | Ticket  | Description                                                          |
|------------|---------|----------------------------------------------------------------------|
| 2018-04-05 | GT-429  | option to turn off exact and graph aligners in grmpy                 |
| 2018-04-05 | GT-428  | upgrade htslib to version 1.8                                        |
| 2018-04-04 | GT-427  | GT-427 multigrmpy to generate graph ID if vc2toparagraph does not provide it|
| 2018-04-03 | GT-415  | GT-415 use parallel grmpy in multigrmpy                              |
| 2018-03-22 | GT-414  | Use filter instead of mapping status to decide whether to graph align after kmer aligner|
| 2018-03-20 | GT-412  | Kmer aligner produces fully clipped nodes at start or end of the alignment|
| 2018-03-20 | GT-407  | Refined edge filter                                                  |
| 2018-03-15 | GT-403  | New threading for grmpy                                              |
| 2018-03-15 | GT-235  | Add documentation for read counting                                  |
| 2018-03-14 | GT-406  | Combine source and sink softclips with adjacent node CIGAR in KmerAligner|
| 2018-03-11 | GT-398  | Minor grmpy fix and documentation update                             |
| 2018-03-07 | GT-394  | graph-level threading to improve efficiency with high-latency file systems|
| 2018-01-03 | GT-393  | Command line option for simulated reads coverage                     |
| 2018-01-03 | GT-392  | Reads sometimes out of order in simulated.bam                        |
| 2018-03-05 | GT-396  | Add command line options to specify location of BAM index separately |
| 2018-02-28 | GT-388  | Improved genotyping on variants with different breakpoint genotypes  |


# Version 1.1

| Date Y-m-d | Ticket  | Description                                                          |
|------------|---------|----------------------------------------------------------------------|
| 2018-02-21 | GT-374  | support for read-level validation                                    |
| 2018-02-19 | GT-379  | configure tool for installation                                      |
| 2018-02-15 | GT-373  | Speedup bam processing by keeping the file open between the graphs   |
| 2018-02-15 | GT-369  | Paragraph can align reads from multiple BAM files                    |
| 2018-02-15 | GT-365  | Add functionality to extract read-supported haplotype paths          |
| 2018-02-09 | GT-360  | Add kmer-based aligner to paragraph                                  |
| 2018-02-09 | GT-310  | Introduce kmer-based readfiltering, improve bad_align filter         |
| 2018-02-07 | GT-294  | Add internal statistics function library, Gaussian fitting           |
| 2018-02-06 | GT-352  | Support in paragraph for running multiple graphs in one go           |
| 2018-02-04 | GT-345  | Refactor CMake sources, split out dependency compile process         |
| 2018-02-01 | GT-306  | Improve build process, allow linking against paragraph as a library  |
| 2018-02-01 | GT-318  | Genotyping parameters from external JSON                             |
| 2018-01-29 | GT-316  | Add graph validation statistics in output JSON                       |
| 2018-01-16 | GT-315  | Refactor genotyping code to allow joint genotyping inside `grmpy`    |
| 2017-12-13 | GT-316  | Add more code for gapless read alignment                             |
| 2017-12-13 | GT-220  | Improve read extraction and increase target region size              |
| 2017-12-11 | GT-304  | Add basic infrastructure for gapless alignments                      |
| 2017-12-07 | GT-305  | Update testing infrastructure + add Docker based tests               |
| 2017-12-06 | GT-311  | Minor build fix for certain versions of g++                          |
| 2017-11-30 | GT-308  | Improved documentation of genotyping model                           |
| 2017-11-28 | GT-292  | Add kmer index class and basic tools for graph indexing              |

# Version 1.0

First release of SV genotyper with support for deletions, insertions, substitutions.
