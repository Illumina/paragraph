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
bin/test_grm
bin/test_blackbox

export PARAGRAPH=/opt/paragraph
export PARAGRAPH_SOURCE=/opt/paragraph-source
export GRMPY_INSTALL=/opt/paragraph

cd ${PARAGRAPH_SOURCE}
python3 ${PARAGRAPH_SOURCE}/src/sh/run_all_python_tests.py