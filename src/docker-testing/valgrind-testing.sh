#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <paragraph-source-directory> <release_tarball_name>"
    exit 1
fi

set -e

if [[ ! -f $2 ]]; then
    echo "Build tarball file does not exist."
    exit 1
fi

realpath() {
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

SOURCE_DIR=$1
RELEASE_DIR=$(dirname $(python ${DIR}/realpath.py $2))
RELEASE_FILE=$(basename $2)

if [[ -z "${RELEASE_DIR// }" ]]; then
    RELEASE_DIR=$(pwd)
fi

docker pull ilmncgrpmi/paragraph-testing-dev:master

docker run --rm  \
    -t \
    -v ${DIR}:/opt/testing-helpers \
    -v ${SOURCE_DIR}:/opt/paragraph-source \
    -v ${RELEASE_DIR}:/release \
    -v $(pwd):/valgrind-test \
     ilmncgrpmi/paragraph-testing-dev:master \
     /bin/bash /opt/testing-helpers/test-helpers/valgrind-test.sh /release/${RELEASE_FILE}