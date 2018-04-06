#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <output_tarball_name>"
    exit 1
fi

set -e

unset HTSLIB_INSTALL_PATH
mkdir /opt/paragraph-build
cd /opt/paragraph-build
cmake /opt/paragraph-source \
    -DCMAKE_C_COMPILER=clang-5.0 \
    -DCMAKE_CXX_COMPILER=clang++-5.0 \
    -DCMAKE_INSTALL_PREFIX=/opt/paragraph && \
    make -j8 && make install

# archive build and install folders
cd /opt
tar czf /workspace/$1-clang-release.tar.gz paragraph/*
tar czf /workspace/$1-clang-buildfolder.tar.gz paragraph-build/*

