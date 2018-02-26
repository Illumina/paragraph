#!/bin/bash
#
# replaces bam header and read alignment information so that simulated reads can be retrieved
# from their position in the whole genome
#
# Author: Roman Petrovski <rpetrovski@illumina.com>


INPUT_BAM=$1
INPUT_FAI=$2

echo -e '@HD\tVN:1.0\tSO:coordinate'

cat ${INPUT_FAI} |awk '{print "@SQ\tSN:"$1"\tLN:"$2}'

samtools view ${INPUT_BAM} |\
awk '{print $3":"(and($2, 16) ? "-" : "+")$4"#"$1"\t"$2"\t"gensub(/.*_(.*)\+.*-.*/, "\\1", $0, $3)"\t"gensub(/.*_.*\+(.*)-.*/, "\\1", $0, $3)"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}'
