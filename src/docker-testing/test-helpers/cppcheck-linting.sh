#!/bin/bash

export PARAGRAPH_SOURCE=/opt/paragraph-source

mkdir -p /opt/paragraph/cppcheck

rm -rf /workspace/cppcheck.xml /workspace/cppcheck-result

cppcheck --enable=information \
         -I${PARAGRAPH_SOURCE}/src/c++/include \
         -I${HTSLIB_INSTALL_PATH}/include \
         --enable=all \
         --suppress=uninitMemberVar \
         --suppress=passedByValue \
         --suppress=ConfigurationNotChecked \
         --suppress=toomanyconfigs \
         --inline-suppr \
         --xml --xml-version=2 \
         --error-exitcode=1 \
         -j 4 \
         ${PARAGRAPH_SOURCE}/src/c++ 2> /workspace/cppcheck.xml 

python ${PARAGRAPH_SOURCE}/src/sh/cppcheck-check.py /workspace/cppcheck.xml

if [[ $? == 0 ]]; then
    echo "No cppcheck issues found."
    mkdir -p /workspace/cppcheck-result
    echo "<html><body>No issues found</body></html>" > /workspace/cppcheck-result/index.html
    exit 0
else
    echo "Cppcheck issues were found, return code was $?. Please inspect cppcheck.xml or cppcheck-result/index.html"
    cppcheck-htmlreport --file=/workspace/cppcheck.xml \
                        --report-dir=/workspace/cppcheck-result >& /dev/null
    exit 1
fi