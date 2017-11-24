#!/usr/bin/env python3
# coding=utf-8

#
# Copyright (c) 2017 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in the LICENSE file in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt
#
#
# June 2017
#
# Convert FASTA MSA to VCF relative to the first sequence
#
# Usage:
#
# For usage instructions run with option --help
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

# noinspection PyUnresolvedReferences
import json
import sys
import gzip
import difflib
from pprint import pformat
from collections import defaultdict

import findgrm  # pylint: disable=unused-import


def read_alignment(filename):
    if filename.endswith("gz"):
        f = gzip.open(filename, "rt")
    else:
        f = open(filename, "rt")
    return json.load(f)


def main():
    data1 = read_alignment(sys.argv[1])
    data2 = read_alignment(sys.argv[2])

    assert "alignments" in data1
    assert "alignments" in data2

    def alignment_pair():
        return {"data1": [], "data2": []}

    alignments_by_id = defaultdict(alignment_pair)

    for aln in data1["alignments"]:
        pos = 0
        chromid = 0
        try:
            chromid = aln["chromId"]
            pos = aln["pos"]
        except:  # pylint: disable=bare-except
            pass
        alid = aln["fragmentId"] + "-mapped-to-%s:%i" % (str(chromid), int(pos))
        alignments_by_id[alid]["data1"].append(aln)

    for aln in data2["alignments"]:
        pos = 0
        chromid = 0
        try:
            chromid = aln["chromId"]
            pos = aln["pos"]
        except:  # pylint: disable=bare-except
            pass
        alid = aln["fragmentId"] + "-mapped-to-%s:%i" % (str(chromid), int(pos))
        alignments_by_id[alid]["data2"].append(aln)

    for k, pair in alignments_by_id.items():
        if len(pair["data1"]) != len(pair["data2"]):
            print("[COUNT] Alignment with different counts: %s / %i != %i" % (k, len(pair["data1"]), len(pair["data2"])))
            continue
        if len(pair["data1"]) != 1:
            print("[ID] Alignment with non-unique ID: %s" % k)
        d1 = pformat(pair["data1"])
        d2 = pformat(pair["data2"])
        if d1 != d2:
            differ = difflib.Differ()
            print("[DIFF] Difference detected for %s:" % k)
            for line in differ.compare(d1.splitlines(), d2.splitlines()):
                print("[DETAILS] " + line)

if __name__ == '__main__':
    main()
