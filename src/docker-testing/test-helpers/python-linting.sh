#!/bin/bash

export PARAGRAPH=/opt/paragraph
export PARAGRAPH_SOURCE=/opt/paragraph-source
export WORKSPACE=/workspace

export PYTHONPATH=${PARAGRAPH_SOURCE}/src/python/lib:${PARAGRAPH_SOURCE}/src/python/bin 

set -e
cd ${WORKSPACE}
python3 -mpylint --rcfile=${PARAGRAPH_SOURCE}/.pylintrc \
    --reports=n \
    ${PARAGRAPH_SOURCE}/src/python/lib \
    ${PARAGRAPH_SOURCE}/src/python/test \
    ${PARAGRAPH_SOURCE}/src/python/bin 

python3 -mpep8 --ignore=E126,E226,E303,E402,E722 --max-line-length=160 \
    ${PARAGRAPH_SOURCE}/src/python