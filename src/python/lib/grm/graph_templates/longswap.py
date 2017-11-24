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


def make_graph(chrom, start, end, ins, flank=150):
    """
    Make a long deletion graph
    :param chrom: chromosome name
    :param start: start coordinate (first deleted base)
    :param end: end coordinate (last deleted base)
    :param ins: inserted sequence
    :param flank: flank length
    :return: paragraph dict
    """
    assert end - start + 1 >= 2 * flank
    target_region_l = "%s:%i-%i" % (chrom,
                                    max(1, start - flank - 1),
                                    start + flank + 1)
    target_region_r = "%s:%i-%i" % (chrom,
                                    max(1, end - flank - 1),
                                    end + flank + 1)

    lf_pos = "%s:%i-%i" % (
        chrom,
        max(1, start - flank - 1),
        start - 1)

    mid_l_pos = "%s:%i-%i" % (
        chrom, start, start + flank - 1)

    mid_r_pos = "%s:%i-%i" % (
        chrom,
        min(1, end - flank),
        min(1, end - 1))

    rf_pos = "%s:%i-%i" % (
        chrom, end + 1, end + flank + 1)

    graph = {
        "sequencenames": ["REF", "DEL", "INS"],
        "target_regions": [target_region_l, target_region_r],
        "nodes": [
            {
                "name": "source",
                "sequence": "NNNNN"
            },
            {
                "name": "LF",
                "reference": lf_pos
            },
            {
                "name": "MID_L",
                "reference": mid_l_pos
            },
            {
                "name": "INS",
                "sequence": ins
            },
            {
                "name": "MID_R",
                "reference": mid_r_pos
            },
            {
                "name": "RF",
                "reference": rf_pos
            },
            {
                "name": "sink",
                "sequence": "NNNNN"
            },
        ],
        "edges": [
            {
                "from": "source",
                "to": "LF",
            },
            {
                "from": "source",
                "to": "MID_R",
            },
            {
                "from": "LF",
                "to": "RF",
                "sequences": ["DEL"]
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
            },
            {
                "from": "LF",
                "to": "MID_L",
                "sequences": ["REF"]
            },
            {
                "from": "MID_R",
                "to": "RF",
                "sequences": ["REF"]
            },
            {
                "from": "MID_R",
                "to": "sink",
            },
            {
                "from": "RF",
                "to": "sink",
            }
        ],
        "paths": [
            {
                "nodes": ["LF", "MID_L"],
                "path_id": "REF|1",
                "sequence": "REF",
                "nucleotide_length": 2 * flank
            },
            {
                "nodes": ["MID_R", "RF"],
                "path_id": "REF|2",
                "sequence": "REF",
                "nucleotide_length": 2 * flank
            },
            {
                "nodes": ["LF", "RF"],
                "path_id": "DEL|1",
                "sequence": "DEL",
                "nucleotide_length": 2 * flank
            },
            {
                "nodes": ["LF", "INS", "RF"],
                "path_id": "INS|1",
                "sequence": "INS",
                "nucleotide_length": 2 * flank + len(ins)
            }
        ]
    }
    return graph
