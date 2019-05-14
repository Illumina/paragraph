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
# Authors:
#
# Felix Schlesinger <fschlesinger@illumina.com>

import logging

from collections import defaultdict

import intervaltree


class GraphContainer:
    """
    Holds a sequence graph in dict format
    """

    def __init__(self, name="VCF Graph"):
        self.name = name
        self.chroms = set()
        self.target_regions = None
        self.sequences = set()
        self.paths = []  # [{path1}, {path2}]
        self.nodes = {}  # nodeName -> node
        self.refNode_starting_at = {}  # (chrom,pos) -> node
        self.refNode_ending_at = {}  # (chrom,pos) -> node
        self.nodes_starting_at = defaultdict(list)  # (chrom,pos) -> [nodes]
        self.edges = {}  # edgeKey -> edge
        self.edges_by_node = defaultdict(list)  # NodeName -> In and out edges of Node

    @staticmethod
    def _edgeKey(n1, n2):
        return n1["name"] + "_" + n2["name"]

    def add_edge(self, nodeFrom, nodeTo, haplos=()):
        logging.debug("Adding edge %s -> %s", nodeFrom['name'], nodeTo['name'])
        ekey = self._edgeKey(nodeFrom, nodeTo)
        if ekey not in self.edges:
            assert nodeFrom["name"] != nodeTo["name"], nodeFrom["name"]
            edge = {
                "from": nodeFrom["name"],
                "to": nodeTo["name"],
                "sequences": set(),
                "name": ekey
            }
            self.edges[ekey] = edge
            self.edges_by_node[nodeFrom["name"]].append(edge)
            self.edges_by_node[nodeTo["name"]].append(edge)
        self.edges[ekey]["sequences"].update(haplos)
        self.sequences.update(haplos)

    def del_edge(self, edge):
        self.edges_by_node[edge["from"]] = [
            e for e in self.edges_by_node[edge["from"]]
            if e["name"] != edge["name"]
        ]
        self.edges_by_node[edge["to"]] = [
            e for e in self.edges_by_node[edge["to"]]
            if e["name"] != edge["name"]
        ]
        del self.edges[edge["name"]]

    def has_edge(self, nodeFrom, nodeTo):
        return self._edgeKey(nodeFrom, nodeTo) in self.edges

    def get_edge(self, nodeNameFrom, nodeNameTo):
        return self.edges[self._edgeKey(self.nodes[nodeNameFrom], self.nodes[nodeNameTo])]

    def inEdges(self, node, haplo=None):
        for e in self.edges_by_node[node["name"]]:
            if e["to"] == node["name"]:
                if haplo is None or haplo in e["sequences"]:
                    yield e

    def outEdges(self, node, haplo=None):
        for e in self.edges_by_node[node["name"]]:
            if e["from"] == node["name"]:
                if haplo is None or haplo in e["sequences"]:
                    yield e

    def altNodes(self, chrom=None):
        for n in self.nodes.values():
            if "reference" not in n:
                if chrom is None or chrom == n.get("chrom"):
                    yield n

    def refNodes(self, chrom=None):
        for n in self.nodes.values():
            if "reference" in n:
                if chrom is None or chrom == n.get("chrom"):
                    yield n

    def add_altNode(self, chrom, start, end, sequence, sequences=(), name=None):
        refSpan = f"{chrom}:{start}-{end}"
        name = name or f"{refSpan}:{sequence}"
        node = {
            "name": name,
            "position": refSpan,
            "sequence": sequence,
            "sequences": set(sequences),
            "chrom": chrom,
            "start": start,
            "end": end,
        }
        self.chroms.add(chrom)
        self.nodes_starting_at[chrom, start].append(node)
        self.nodes[name] = node
        self.sequences.update(sequences)
        logging.debug("Created ALT node: %s", str(node))
        return node

    def add_refNode(self, chrom, start, end, sequences=(), name=None):
        refSpan = f"{chrom}:{start}-{end}"
        node = {
            "name": name or f"ref-{refSpan}",
            "reference": refSpan,
            "sequences": set(sequences),
            "chrom": chrom,
            "start": start,
            "end": end,
        }
        self.chroms.add(chrom)
        self.refNode_starting_at[chrom, start] = node
        self.refNode_ending_at[chrom, end] = node
        self.nodes_starting_at[chrom, start].append(node)
        self.nodes[node["name"]] = node
        self.sequences.update(sequences)
        logging.debug("Created REF node: %s", str(node))
        return node

    def del_node(self, node):
        """Delete node and all its edges from graph"""
        for e in self.edges_by_node[node["name"]]:
            self.del_edge(e)
        self.nodes_starting_at[node["start"]] = [
            n for n in self.nodes_starting_at[node["start"]]
            if n["name"] != node["name"]
        ]
        del self.nodes[node["name"]]

    def nodes_by_haplo(self, haplo):
        ns = [n for n in self.nodes.values()
              if haplo in n["sequences"]]
        ns.sort(key=lambda n: (n["start"], n["end"]))
        return ns

    def check(self):
        for e in self.edges.values():
            assert e["from"] in self.nodes
            assert e["to"] in self.nodes
            assert e["from"] != e["to"]
            assert all(s in self.sequences for s in e.get("sequence", ()))
        for p in self.paths:
            assert p["sequence"] in self.sequences
            for n in p["nodes"]:
                assert n in self.nodes

    def topological_sort(self):
        """
        Get nodes and edges in topological sort order
        :return: sorted_nodes, sorted_edges
        """
        ordered_nodes = []
        for n in self.nodes.values():
            if "mark" in n:
                del n["mark"]

        def visit(node):
            node["mark"] = "t"
            nnodes = [self.nodes[e["to"]] for e in self.outEdges(node)]
            for x in sorted(nnodes, key=lambda n: n["name"]):
                if "mark" not in x:
                    visit(x)
                elif x["mark"] == "t":
                    raise Exception("Graph has a cycle at %s -> %s" % (node["name"], x["name"]))
            node["mark"] = "p"
            ordered_nodes.insert(0, node)

        for node in self.nodes.values():
            if "mark" not in node:
                visit(node)
        nodeOrder = {n["name"]: i for (i, n) in enumerate(ordered_nodes)}
        all_edges = sorted(self.edges.values(),
                           key=lambda e: (nodeOrder[e["from"]], nodeOrder[e["to"]]))
        return ordered_nodes, all_edges

    def get_reference_regions(self):
        """
        Genome regions covered by reference nodes
        """
        for chrom in self.chroms:
            tree = intervaltree.IntervalTree()
            for n in self.refNodes(chrom):
                tree.addi(n["start"], n["end"] + 1)
            tree.merge_overlaps()
            for iv in tree.items():
                yield f"{chrom}:{iv.begin}-{iv.end-1}"

    def json_dict(self):
        """
        Dictionary with graph in paragraph json representation
        """
        attribsToDel = ["mark", "vcfId", "chrom", "start", "end"]
        sNodes, sEdges = self.topological_sort()
        nodes = []
        for n in sNodes:
            n = n.copy()
            for a in attribsToDel + ["sequences"]:
                if a in n:
                    del n[a]
            nodes.append(n)
        edges = []
        for e in sEdges:
            e = e.copy()
            for a in attribsToDel:
                if a in e:
                    del e[a]
            if "sequences" in e:
                e["sequences"] = list(sorted(e["sequences"]))
                if not e["sequences"]:
                    del e["sequences"]
            edges.append(e)

        return {
            "nodes": nodes,
            "edges": edges,
            "paths": self.paths,
            "target_regions": list(sorted(self.target_regions)),
            "sequencenames": list(sorted(self.sequences)),
            "model_name": self.name
        }
