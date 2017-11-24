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


def make_graph(chrom, start, end, flank=150):
    """
    Make a small deletion graph
    :param chrom: chromosome name
    :param start: start coordinate (first deleted base)
    :param end: end coordinate (last deleted base)
    :param flank: flank length
    :return: paragraph dict
    """
    target_region = "%s:%i-%i" % (chrom,
                                  max(1, start - flank - 1),
                                  end + flank + 1)

    lf_pos = "%s:%i-%i" % (
        chrom,
        max(1, start - flank - 1),
        start - 1)

    mid_pos = "%s:%i-%i" % (
        chrom, start, end)

    rf_pos = "%s:%i-%i" % (
        chrom, end + 1, end + flank + 1)

    graph = {
        "sequencenames": ["REF", "DEL"],
        "target_regions": [target_region],
        "nodes": [
            {
                "name": "LF",
                "reference": lf_pos
            },
            {
                "name": "MID",
                "reference": mid_pos
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
                "sequences": ["DEL"]
            },
            {
                "from": "LF",
                "to": "MID",
                "sequences": ["REF"]
            },
            {
                "from": "MID",
                "to": "RF",
                "sequences": ["REF"]
            }
        ],
        "paths": [
            {
                "nodes": ["LF", "MID", "RF"],
                "path_id": "REF|1",
                "sequence": "REF",
                "nucleotide_length": end - start + 1 + 2 * flank
            },
            {
                "nodes": ["LF", "RF"],
                "path_id": "DEL|1",
                "sequence": "DEL",
                "nucleotide_length": 2 * flank
            }
        ]
    }
    return graph
