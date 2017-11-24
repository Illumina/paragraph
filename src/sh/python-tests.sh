#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [[ -f /illumina/development/grmpy/py3-ve/bin/activate ]]; then
    . /illumina/development/grmpy/py3-ve/bin/activate
elif [[ -f /bioinfoSD/users/pkrusche/py3/bin/activate ]]; then
    . /bioinfoSD/users/pkrusche/py3/bin/activate
fi

python3 ${DIR}/run_all_python_tests.py