#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <paragraph-source-directory> <clang_build_folder_tarball>"
    exit 1
fi

docker pull ilmncgrpmi/paragraph-testing-dev:clang-master

realpath() {
    [[ $1 = /* ]] && echo "$1" || echo "$PWD/${1#./}"
}

SOURCE_DIR=$1
BUILD_FOLDER_DIR=$(dirname $(python ${DIR}/realpath.py $2))
BUILD_FOLDER_FILE=$(basename $2)

docker run --rm  \
    -t \
    -v $DIR:/opt/testing-helpers \
    -v $1:/opt/paragraph-source \
    -v `pwd`:/workspace \
    -v ${BUILD_FOLDER_DIR}:/build-folder \
     ilmncgrpmi/paragraph-testing-dev:clang-master \
     /bin/bash /opt/testing-helpers/test-helpers/clang-linting.sh /build-folder/${BUILD_FOLDER_FILE}

