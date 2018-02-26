#!/bin/bash

export PARAGRAPH_SOURCE=/opt/paragraph-source
export CLANG_FORMAT=clang-format-5.0

set +e
DIFFS=0
for x in $(find ${PARAGRAPH_SOURCE}/src/c++ -iname *.h -o -iname *.hh -o -iname *.cpp); do
    ${CLANG_FORMAT} -style=file $x | diff $x -
    if [[ $? != 0 ]]; then
        echo "Differences found in $x"
        DIFFS=1
    fi
done

if [[ $DIFFS == 0 ]]; then
    echo "No formatting issues found."
    exit 0
else
    echo "Some formatting issues were found."
    exit 1
fi
