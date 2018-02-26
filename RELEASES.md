# Paragraph Release Notes / Change Log

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
