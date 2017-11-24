#!/usr/bin/env bash

set -e

WORKSPACE=$(pwd)

. ${WORKSPACE}/src/sh/illumina-setup.sh

cd ${WORKSPACE}/install

# google's death tests break valgrind, so we exclude them here
valgrind --leak-check=full --xml=yes \
		 --xml-file=${WORKSPACE}/valgrind1.xml \
         --suppressions=${WORKSPACE}/src/sh/valgrind-suppressions.supp \
           ${WORKSPACE}/install/bin/test_grm --gtest_filter=-*DeathTest

python ${WORKSPACE}/src/sh/valgrind-check.py ${WORKSPACE}/valgrind1.xml

export HG19=/illumina/sync/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa
export HG38=/illumina/sync/igenomes/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa

valgrind --leak-check=full --xml=yes \
		 --xml-file=${WORKSPACE}/valgrind2.xml \
         --suppressions=${WORKSPACE}/src/sh/valgrind-suppressions.supp \
           ${WORKSPACE}/install/bin/test_blackbox

python ${WORKSPACE}/src/sh/valgrind-check.py ${WORKSPACE}/valgrind2.xml

valgrind --leak-check=full --xml=yes \
		 --xml-file=${WORKSPACE}/valgrind3.xml \
         --suppressions=${WORKSPACE}/src/sh/valgrind-suppressions.supp \
           ${WORKSPACE}/install/bin/grmpy \
           -r  ${WORKSPACE}/share/test-data/paragraph/long-del/chr4_graph_typing.fa \
           -p  ${WORKSPACE}/share/test-data/genotyping_test/chr4_graph_typing.2sample.json \
           -m ${WORKSPACE}/share/test-data/genotyping_test/chr4_graph_typing.manifest \
           -o t${WORKSPACE}/vg_test.json

python ${WORKSPACE}/src/sh/valgrind-check.py ${WORKSPACE}/valgrind3.xml
