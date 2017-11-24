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
# August 2017
#
# Split long reference nodes.
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

from copy import copy
from grm.helpers import parse_region
from grm.vcfgraph.topological_sort import topological_sort


def split_long_reference_nodes(all_nodes, all_edges,
                               max_len=300, padding_len=150):
    """
    Split long reference nodes.
    :param all_nodes: list of nodes
    :param all_edges: list of edges
    :param max_len: max length of reference node with no sequences
    :param padding_len: length of sequence to keep
    :return: nodes, edges
    """
    assert max_len >= 2 * padding_len
    new_nodes = []
    new_edges = copy(all_edges)
    for i, n in enumerate(all_nodes):
        # only break ref-nodes that don't specify a sequence
        if "reference" not in n or ("sequences" in n and (n["sequences"] != ["REF"] and n["sequences"])):
            new_nodes.append(n)
            continue
        chrom, start_pos, end_pos = parse_region(n["reference"])
        if end_pos - start_pos + 1 <= max_len:
            new_nodes.append(n)
            continue
        node1 = {
            "name": "ref-%s:%i-%i-%i" % (chrom, start_pos, start_pos + padding_len - 1, i),
            "reference": "%s:%i-%i" % (chrom, start_pos, start_pos + padding_len - 1),
            "start": start_pos,
            "end": start_pos + padding_len - 1
        }
        node2 = {
            "name": "ref-%s:%i-%i-%i" % (chrom, end_pos - padding_len + 1, end_pos, i),
            "reference": "%s:%i-%i" % (chrom, end_pos - padding_len + 1, end_pos),
            "start": end_pos - padding_len + 1,
            "end": end_pos
        }
        new_nodes += [node1, node2]

        relevant_edges = [e for e in new_edges if e["from"] == n["name"] or e["to"] == n["name"]]

        for e in relevant_edges:
            if e["from"] == n["name"]:
                e["from"] = node2["name"]
            elif e["to"] == n["name"]:
                e["to"] = node1["name"]

    try:
        if new_nodes[0]["sequence"] == "N" * new_nodes[0]["sequence"]:
            source_name = new_nodes[0]["name"]
        else:
            raise ValueError()
    except:  # pylint: disable=bare-except
        source_name = "source"

    try:
        if new_nodes[-1]["sequence"] == "N" * new_nodes[-1]["sequence"]:
            sink_name = new_nodes[-1]["name"]
        else:
            raise ValueError()
    except:  # pylint: disable=bare-except
        sink_name = "sink"

    new_nodes, new_edges = topological_sort(new_nodes, new_edges, source_name, sink_name)

    return new_nodes, new_edges
