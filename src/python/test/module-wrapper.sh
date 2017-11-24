#!/usr/bin/env bash

if [[ -d /illumina ]]; then
    module purge
fi

$@