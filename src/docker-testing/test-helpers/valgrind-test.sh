#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <build_tarball_name>"
    exit 1
fi

if [[ ! -f $1 ]]; then
    echo "Build tarball file $1 does not exist in Docker image."
    exit 1
fi

set -e
cd /opt
tar xf $1
cd paragraph

mkdir -p /valgrind-test

export PARAGRAPH=/opt/paragraph
export PARAGRAPH_SOURCE=/opt/paragraph-source
export WORKSPACE=/valgrind-test

cd $WORKSPACE

valgrind --leak-check=full --xml=yes \
         --xml-file=${WORKSPACE}/valgrind1.xml \
         --suppressions=${PARAGRAPH_SOURCE}/src/sh/valgrind-suppressions.supp \
           ${PARAGRAPH}/bin/test_grm \
         --gtest_filter=-*Performance*:IntervalTree*:-*Performance*:*Random

python ${PARAGRAPH_SOURCE}/src/sh/valgrind-check.py ${WORKSPACE}/valgrind1.xml

valgrind --leak-check=full --xml=yes \
		 --xml-file=${WORKSPACE}/valgrind2.xml \
         --suppressions=${PARAGRAPH_SOURCE}/src/sh/valgrind-suppressions.supp \
           ${PARAGRAPH}/bin/test_blackbox

python ${PARAGRAPH_SOURCE}/src/sh/valgrind-check.py ${WORKSPACE}/valgrind2.xml

valgrind --leak-check=full --xml=yes \
		 --xml-file=${WORKSPACE}/valgrind3.xml \
         --suppressions=${PARAGRAPH_SOURCE}/src/sh/valgrind-suppressions.supp \
           ${PARAGRAPH}/bin/grmpy \
           -r  ${PARAGRAPH}/share/test-data/paragraph/long-del/chr4_graph_typing.fa \
           -g  ${PARAGRAPH}/share/test-data/paragraph/long-del/chr4_graph_typing.2sample.json \
           -m ${PARAGRAPH}/share/test-data/paragraph/long-del/chr4_graph_typing.manifest \
           -o ${WORKSPACE}/vg_test.json

python ${PARAGRAPH_SOURCE}/src/sh/valgrind-check.py ${WORKSPACE}/valgrind3.xml

