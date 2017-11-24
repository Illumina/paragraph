#!/bin/bash

set -e

# Find python
PYTHON=python
WD=$(pwd)
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "Checking compiled dependencies..."

if [ "$1" == "clean" ]; then
    rm -rf ${WD}/protobuf-install ${WD}/htslib
    exit 0
fi

# wget https://github.com/google/protobuf/releases/download/v3.0.0-beta-1/protobuf-cpp-3.0.0-beta-1.tar.gz

if [[ -z ${PROTOBUF_INSTALL_PATH} || ! -d ${PROTOBUF_INSTALL_PATH} ]]; then
    if [ ! -d ${WD}/protobuf-install ] ;
    then
        echo "[Building protobuf]"
        tar xzf ${DIR}/protobuf-cpp-3.3.0.tar.gz
        cd protobuf-3.3.0
        ./autogen.sh
        PKG_CONFIG_PATH=${WD}/protobuf-install/lib/pkgconfig ./configure --prefix=$WD/protobuf-install
        make -j4
        make -j4 install
    else
        echo "Protobuf already built. To rebuild, delete ${WD}/protobuf-install"
    fi
else
    echo "Using pre-built protobuf from ${PROTOBUF_INSTALL_PATH}"
fi

# Build htslib

if [ ! -d ${WD}/htslib ] ;
then
    echo "[Building HTSlib]"
    mkdir ${WD}/htslib
    cd ${WD}/htslib
    tar xzf ${DIR}/htslib.tar.gz
    make -j4
else
    echo "HTSlib is already built. To rebuild, delete ${WD}/htslib"
fi

# Build spdlog

if [ ! -d ${WD}/spdlog ] ;
then
    echo "[Extracting spdlog]"
    mkdir -p ${WD}/spdlog
    cd ${WD}/spdlog
    tar xzf ${DIR}/spdlog.tar.gz
else
    echo "Spdlog is already extracted. To rebuild, delete ${WD}/spdlog"
fi

# Build Google Test/Mock

if [ ! -d ${WD}/googletest-release-1.8.0 ] ;
then
    echo "[Extracting googletest]"
    cd ${WD}
    tar xzf ${DIR}/googletest-release-1.8.0.tar.gz
else
    echo "Google Test is already extracted. To built, delete ${WD}/googletest-release-1.8.0"
fi