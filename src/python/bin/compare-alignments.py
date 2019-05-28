#!/usr/bin/env python3

# coding=utf-8
#
# Copyright (c) 2016-2019 Illumina, Inc.
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# You may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied
# See the License for the specific language governing permissions and limitations
#
import sys
import difflib
from pprint import pformat
from collections import defaultdict

import findgrm  # pylint: disable=unused-import
from grm.helpers import load_json


def main():
    data1 = load_json(sys.argv[1])
    data2 = load_json(sys.argv[2])
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
