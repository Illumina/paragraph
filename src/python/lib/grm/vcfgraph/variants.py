# coding=utf-8
#
# Copyright (c) 2018 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in the LICENSE file in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt
#
# Author:
#
# Felix Schlesinger <fschlesinger@illumina.com>

import sys
import logging
from collections import namedtuple
import traceback

from grm.vcfgraph.graphContainer import GraphContainer


Variant = namedtuple("Variant", ["node", "start", "end", "alt"])


def add_variants(graph: GraphContainer, variants):
    """
    Add variants to a graph as new nodes and edges
    :param graph: The graph object to modify
    :param variants: 'variants' entry from the paragraph output JSON
    """
    for (node, nodeVarDicts) in variants.items():
        nodeVars = []
        for varDict in nodeVarDicts:
            # older / protobuf based outputs don't have these when empty
            start_pos = varDict.get("start", 0)
            end_pos = varDict.get("end", 0)
            alt = varDict.get("alt", "")
            nodeVars.append(Variant(node, start_pos, end_pos, alt))
        try:
            add_variants_node(graph, graph.nodes[node], nodeVars)
        except Exception:  # pylint: disable=broad-except
            traceback.print_exc(file=sys.stderr)
            print("Skipping variant records on node {}: {}.".format(node, nodeVars), file=sys.stderr)


def split_node(graph: GraphContainer, node, breakpoints):
    """
    Split a node at a set of breakpoints and link new (sub-)nodes
    Used to link to new variant nodes later
    Modifies graph and deletes node after splitting
    :returns Created sub-nodes
    """
    if not breakpoints:
        return node
    breakpoints = sorted(set(breakpoints))
    logging.debug("Splitting %s at %s ", node['name'], str(breakpoints))
    nodes = []
    lEnd = 0
    for p in breakpoints:
        assert 0 <= p <= node["end"] - node["start"] + 1
        nStart = node["start"] + lEnd
        nEnd = node["start"] + p - 1
        if "reference" in node:
            nodes.append(graph.add_refNode(node["chrom"], nStart, nEnd, node["sequences"]))
        else:
            seq = node["sequence"][lEnd:p]
            nodes.append(graph.add_altNode(node["chrom"], nStart, nEnd, seq, node["sequences"]))
        lEnd = p
    # Add last node
    lStart = node["start"] + breakpoints[-1]
    if "reference" in node:
        nodes.append(graph.add_refNode(node["chrom"], lStart, node["end"], node["sequences"]))
    else:
        seq = node["sequence"][breakpoints[-1]:]
        nodes.append(graph.add_altNode(node["chrom"], lStart, node["end"], seq, node["sequences"]))
    # Connect nodes
    for e in graph.inEdges(node):
        graph.add_edge(graph.nodes[e["from"]], nodes[0], e["sequences"])
    for e in graph.outEdges(node):
        graph.add_edge(nodes[-1], graph.nodes[e["to"]], e["sequences"])
    for (n1, n2) in zip(nodes[:-1], nodes[1:]):
        graph.add_edge(n1, n2)
    # Delete original node, unless identical to new node (no split)
    if node['name'] not in [n['name'] for n in nodes]:
        graph.del_node(node)
    return nodes


def add_variants_node(graph: GraphContainer, node, variants):
    """ Add variants to one node in the graph """
    bps = []
    for var in variants:
        if var.start > var.end + 1:
            raise Exception("Variant start({}) > variane end({})!".format(var.start, var.end))
        if var.start == var.end + 1 and not var.alt:
            raise Exception("Variant start({}) == end but no insertion sequence is specified.".format(var.start))
        bps.extend((var.start, var.end+1))
    nodes = split_node(graph, node, bps)
    nodesEnding = {node["end"]: node for node in nodes[:-1]}
    nodesStarting = {node["start"]: node for node in nodes}
    for var in variants:
        vStart = node["start"] + var.start
        vEnd = node["start"] + var.end
        alt = graph.add_altNode(node["chrom"], vStart, vEnd, var.alt)
        graph.add_edge(nodesEnding[vStart - 1], alt)
        graph.add_edge(alt, nodesStarting[vEnd + 1])
