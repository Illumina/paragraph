#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <clang_build_folder_tarball>"
    exit 1
fi

# extract build folder
set -e
cd /opt
tar xf $1
cd paragraph-build

export PARAGRAPH=/opt/paragraph-build
export PARAGRAPH_SOURCE=/opt/paragraph-source

export CLANG_TIDY=clang-tidy-5.0

ISSUES=0
for x in $(find ${PARAGRAPH_SOURCE}/src/c++ -iname *.cpp); do
    echo "Clang-tidy check for $x..."
    ${CLANG_TIDY} -p ${PARAGRAPH} $x
    if [[ $? != 0 ]]; then
        echo "Issues found in $x"
        ISSUES=1
    fi
done

if [[ $ISSUES == 0 ]]; then
    echo "No linting issues found."
    exit 0
else
    echo "Some linting issues were found."
    exit 1
fi