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
# October 2017
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#


def make_graph(chrom, start, ins, flank=150):
    """
    Make a small deletion graph
    :param chrom: chromosome name
    :param start: start coordinate (last base before insertion)
    :param ins: inserted sequence
    :param flank: flank length
    :return: paragraph dict
    """
    target_region = "%s:%i-%i" % (chrom,
                                  max(1, start - flank - 1),
                                  start + flank + 1)

    lf_pos = "%s:%i-%i" % (
        chrom,
        max(1, start - flank - 1),
        max(1, start - 1))

    rf_pos = "%s:%i-%i" % (
        chrom, start + 1, start + flank + 1)

    graph = {
        "sequencenames": ["REF", "INS"],
        "target_regions": [target_region],
        "nodes": [
            {
                "name": "LF",
                "reference": lf_pos
            },
            {
                "name": "INS",
                "sequence": ins
            },
            {
                "name": "RF",
                "reference": rf_pos
            },
        ],
        "edges": [
            {
                "from": "LF",
                "to": "RF",
                "sequences": ["REF"]
            },
            {
                "from": "LF",
                "to": "INS",
                "sequences": ["INS"]
            },
            {
                "from": "INS",
                "to": "RF",
                "sequences": ["INS"]
            }
        ],
        "paths": [
            {
                "nodes": ["LF", "INS", "RF"],
                "path_id": "INS|1",
                "sequence": "INS",
                "nucleotide_length": len(ins) + 2 * flank
            },
            {
                "nodes": ["LF", "RF"],
                "path_id": "REF|1",
                "sequence": "REF",
                "nucleotide_length": 2 * flank
            }
        ]
    }
    return graph
