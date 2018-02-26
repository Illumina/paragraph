#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <paragraph-source-directory> <output_tarball_name>"
    echo "       The output tarball name should be a filename prefix within the current directory."
    exit 1
fi

docker pull ilmncgrpmi/paragraph-testing-dev:clang-master

docker run --rm  \
    -t \
    -v $DIR:/opt/testing-helpers \
    -v $1:/opt/paragraph-source \
    -v `pwd`:/workspace \
     ilmncgrpmi/paragraph-testing-dev:clang-master \
     /bin/bash /opt/testing-helpers/test-helpers/build_clang.sh $2


