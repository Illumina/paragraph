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
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
# Felix Schlesinger <fschlesinger@illumina.com>

import logging

from grm.vcfgraph.graphContainer import GraphContainer
from grm.helpers import parse_region


def add_source_sink(graph: GraphContainer,
                    source_name="source",
                    sink_name="sink"):
    """
    add source and sink if necessary and link to exising nodes without in/out edges
    :param graph: graph to work on
    :param source_name: name of source node
    :param sink_name: name of sink node
    """
    if source_name not in graph.nodes:
        graph.nodes[source_name] = {
            "name": source_name,
            "sequence": "N" * 10
        }
    if sink_name not in graph.nodes:
        graph.nodes[sink_name] = {
            "name": sink_name,
            "sequence": "N" * 10
        }
    # Link nodes without incoming/outgoing edges to source/sink
    for node in graph.nodes.values():
        if node["name"] in [source_name, sink_name]:
            continue
        if not any(graph.inEdges(node)):
            logging.info("Linking %s from source", node['name'])
            graph.add_edge(graph.nodes[source_name], node)
        if not any(graph.outEdges(node)):
            logging.info("Linking %s to sink", node['name'])
            graph.add_edge(node, graph.nodes[sink_name])


def split_ref_nodes(graph: GraphContainer, max_len=300, padding_len=150):
    """
    Split long reference nodes
    :param graph: graph to work on
    :param max_len: max length of reference node with no sequences
    :param padding_len: length of sequence to keep
    """
    assert max_len >= 2 * padding_len
    for node in list(graph.refNodes()):
        if node["end"] - node["start"] + 1 <= max_len:
            continue
        logging.info("Splitting long REF node: %s", node['name'])
        firstEnd = node["start"] + padding_len - 1
        n1 = graph.add_refNode(node["chrom"], node["start"], firstEnd, node["sequences"])
        sndStart = node["end"] - padding_len + 1
        n2 = graph.add_refNode(node["chrom"], sndStart, node["end"], node["sequences"])

        for e in list(graph.inEdges(node)):
            graph.add_edge(graph.nodes[e["from"]], n1, e["sequences"])
        for e in list(graph.outEdges(node)):
            graph.add_edge(n2, graph.nodes[e["to"]], e["sequences"])
        graph.del_node(node)


def split_alt_nodes(graph: GraphContainer, max_len=300, padding_len=150):
    """
    Split long alternate nodes
    :param graph: graph to work on
    :param max_len: max length of reference node with no sequences
    :param padding_len: length of sequence to keep
    """
    assert max_len >= 2 * padding_len
    for node in list(graph.altNodes()):
        if len(node["sequence"]) <= max_len:
            continue
        logging.info("Splitting long ALT node: %s", node['name'])

        n1 = graph.add_altNode(node["chrom"], node["start"], node["end"],
                               node["sequence"][:padding_len], node["sequences"])
        n2 = graph.add_altNode(node["chrom"], node["start"], node["end"],
                               node["sequence"][-padding_len:], node["sequences"])

        for e in list(graph.inEdges(node)):
            graph.add_edge(graph.nodes[e["from"]], n1, e["sequences"])
        for e in list(graph.outEdges(node)):
            graph.add_edge(n2, graph.nodes[e["to"]], e["sequences"])
        graph.del_node(node)


def remove_empty_nodes(graph: GraphContainer):
    """
    Remove nodes without sequence (from deletions / skipped insertions or split ref nodes)
    Merge in & out edge pairs to keep connections
    """
    for node in list(graph.nodes.values()):
        if (("reference" in node and node["start"] <= node["end"]) or
                node.get("sequence", "") != ""):
            continue
        logging.info("Removing empty node %s", node['name'])
        inSeqs = [s for e in graph.inEdges(node) for s in e["sequences"]]
        outSeqs = [s for e in graph.outEdges(node) for s in e["sequences"]]
        for e1 in list(graph.inEdges(node)):
            for e2 in list(graph.outEdges(node)):
                # Label the new edges with sequence labels either observed
                # on both merged in- and out-edge or on an in (out) -edge only
                # if the label is undetermined goung out (in)
                haplos = e1["sequences"].intersection(e2["sequences"]).union(
                    e1["sequences"].difference(outSeqs).union(
                        e2["sequences"].difference(inSeqs)))
                graph.add_edge(graph.nodes[e1["from"]], graph.nodes[e2["to"]], haplos)
        graph.del_node(node)


def combine_nodes(graph: GraphContainer):
    """
    Combine adjacent nodes with the same sequence labels
    """
    for n1 in list(graph.nodes.values()):
        if len(list(graph.outEdges(n1))) != 1:
            continue  # Pair of nodes with no other in/out edges
        n2 = graph.nodes[next(graph.outEdges(n1))["to"]]
        if len(list(graph.inEdges(n2))) != 1:
            continue
        if not (n1["chrom"] == n2["chrom"] and n1["end"] + 1 == n2["start"]):
            continue  # nodes must be adjacent
        haplos = n1["sequences"]
        if n2["sequences"] != haplos:
            continue  # only collapse nodes with same haplotypes
        if "reference" in n1:
            if "reference" not in n2:
                continue  # nodes must be of same type
            node = graph.add_refNode(n1["chrom"], n1["start"], n2["end"], haplos)
        else:
            if "reference" in n2:
                continue  # nodes must be of same type
            node = graph.add_altNode(
                n1["chrom"], n1["start"], n2["end"], n1["sequence"] + n2["sequence"], haplos)
        logging.info("Combinding %s and %s", n1['name'], n2['name'])
        for e in list(graph.inEdges(n1)):
            graph.add_edge(graph.nodes[e["from"]], node, e["sequences"])
        for e in list(graph.outEdges(n2)):
            graph.add_edge(node, graph.nodes[e["to"]], e["sequences"])
        graph.del_node(n1)
        graph.del_node(n2)


def remove_redundant_edge_labels(graph: GraphContainer):
    """
    Remove edge sequence labels that don't add information about haplotypes.
    """
    # If all out-edges of a node have the same edge label
    # and the node is already labeled on an in-edge
    # and those edges to not connect to a node with that label
    for node in graph.nodes.values():
        for haplo in node["sequences"]:
            for e in graph.outEdges(node):
                if haplo not in e["sequences"]:
                    break
                if haplo in graph.nodes[e["to"]]["sequences"]:
                    break
            else:
                for e in graph.outEdges(node):
                    e["sequences"].remove(haplo)


def get_path(graph: GraphContainer, sequence):
    """
    Return paths (list of nodes) covering all edges for one sequence (haplotype)
    :param graph: graph to work on
    :param sequence: Haplotype to cover
    """
    nodes, edges = graph.topological_sort()
    for e in edges:
        if "mark" in e:
            del e["mark"]
    paths = []

    def visit(edge, curPath):
        node = graph.nodes[edge["to"]]
        curPath = curPath + [node["name"]]
        edge["mark"] = True
        paths = []
        for e in graph.outEdges(node, sequence):
            if "mark" not in e:
                paths.extend(visit(e, curPath))
        if not paths:
            paths = [curPath]
        return paths

    paths = []
    for node in nodes:
        for edge in graph.outEdges(node, sequence):
            if "mark" not in edge:
                paths += visit(edge, [node["name"]])
    return paths


def ref_paths(graph):
    for fNode in graph.refNodes():
        for edge in graph.outEdges(fNode):
            nNode = graph.nodes[edge["to"]]
            if "reference" in nNode and fNode["end"] + 1 == nNode["start"]:
                graph.add_edge(fNode, nNode, ["REF"])
    res = []
    for path in get_path(graph, "REF"):
        res.append({
            "nodes": path,
            "path_id": f"REF|{len(res)+1}",
            "sequence": "REF"
        })
    return res


def add_ref_path(graph: GraphContainer):
    """
    Add 'REF' haplotype to all REF edges.
    Create paths covering all 'REF' edges
    """
    for path in ref_paths(graph):
        graph.paths.append(path)


def add_alt_paths(graph: GraphContainer):
    res = []
    refPaths = ref_paths(graph)
    for path in get_path(graph, None):
        if path[0] == "source":
            path = path[1:]
        if path[-1] == "sink":
            path = path[:-1]
        if path not in [p["nodes"] for p in refPaths]:
            res.append({
                "nodes": path,
                "path_id": f"ALT|{len(res)+1}",
                "sequence": "ALT",
            })
            graph.sequences.add("ALT")
    graph.paths += res


def add_graph(graph1: GraphContainer, graph2: GraphContainer):
    """
    Add all nodes, edges and paths from graph2 to graph1 (inplace)
    """
    for node in graph2.refNodes():
        graph1.add_refNode(
            node["chrom"], node["start"], node["end"], node["sequences"])
    for node in graph2.altNodes():
        graph1.add_altNode(
            node["chrom"], node["start"], node["end"], node["sequence"], node["sequences"])
    for edge in graph2.edges.values():
        graph1.add_edge(graph1.nodes[edge["from"]], graph1.nodes[edge["to"]], edge["sequences"])
    graph1.paths += graph2.paths


def load_json(json) -> GraphContainer:
    """
    Construct graph object from JSON representation
    :param json: Dictionary of JSON file contents
    """
    graph = GraphContainer()
    for node in json["nodes"]:
        seqs = node.get("sequences", ())
        if "reference" in node:
            chrom, start, end = parse_region(node["reference"])
            graph.add_refNode(chrom, start, end, seqs, node["name"])
        elif "position" in node:
            chrom, start, end = parse_region(node["position"])
            graph.add_altNode(chrom, start, end, node["sequence"], seqs, node["name"])
        else:
            graph.nodes[node["name"]] = node
    for edge in json["edges"]:
        seqs = edge.get("sequences", ())
        graph.add_edge(graph.nodes[edge["from"]], graph.nodes[edge["to"]], seqs)
    graph.name = json["model_name"]
    graph.paths = json.get("paths", [])
    graph.target_regions = json.get("target_regions", [])
    graph.check()
    return graph
