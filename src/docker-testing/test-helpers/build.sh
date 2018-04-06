#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <output_tarball_name>"
    exit 1
fi

set -e

unset HTSLIB_INSTALL_PATH
mkdir /opt/paragraph-build
cd /opt/paragraph-build
cmake /opt/paragraph-source -DCMAKE_INSTALL_PREFIX=/opt/paragraph && make -j4 && make install

# archive build and install folders
cd /opt
tar czf /workspace/$1-gcc-release.tar.gz paragraph/*
tar czf /workspace/$1-gcc-buildfolder.tar.gz paragraph-build/*

