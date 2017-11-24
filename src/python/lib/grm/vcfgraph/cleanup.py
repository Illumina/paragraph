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
# Cleanup nodes and edges
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

from copy import copy
from collections import defaultdict

from grm.helpers import parse_region


def contract_paths(all_nodes, all_edges):
    """
    Contract adjacent nodes with the same sequence labels into one
    """
    fadj = defaultdict(list)
    radj = defaultdict(list)

    for e in all_edges:
        fadj[e["from"]].append(e["to"])
        radj[e["to"]].append(e["from"])

    nodes_by_name = {n["name"]: n for n in all_nodes}
    edges_by_name = {e["from"] + "_" + e["to"]: e for e in all_edges}

    # preserve order
    for i, n in enumerate(all_nodes):
        n["topo_order"] = i

    nodes_to_check = list([n for n in all_nodes
                           if len(fadj[n["name"]]) == 1 and len(radj[fadj[n["name"]][0]]) == 1])

    while nodes_to_check:
        n1 = nodes_to_check.pop(0)["name"]
        n2 = fadj[n1][0]

        # is n1 the only predecessor?
        if len(radj[n2]) != 1:
            continue

        n1i = nodes_by_name[n1]
        n2i = nodes_by_name[n2]

        # nodes must be of same type
        if "sequence" in n1i and "reference" in n2i or "reference" in n1i and "sequence" in n2i:
            continue

        # reference nodes must be adjacent
        if "reference" in n1i:
            chrom1, start1, end1 = parse_region(n1i["reference"])
            chrom2, start2, end2 = parse_region(n2i["reference"])
            if chrom1 != chrom2 or end1 != start2 - 1:
                continue
        else:
            chrom1, start1, end1 = None, None, None
            chrom2, start2, end2 = None, None, None

        if "sequences" in n1i and "sequences" in n2i:
            if set(n1i["sequences"]) != set(n2i["sequences"]):
                continue
        elif "sequences" in n1i or "sequences" in n2i:
            continue

        # does the edge in between add any disambig sequences?
        e = edges_by_name[n1 + "_" + n2]
        if "sequences" in n1i and "sequences" in n2i:
            if "sequences" in e and set(e["sequences"]) != set(n1i["sequences"]):
                continue

        new_node = copy(n1i)
        if "reference" in n1i:
            new_node["reference"] = "%s:%i-%i" % (chrom1, start1, end2)
            new_node["name"] = "ref-" + new_node["reference"]
        else:  # sequence nodes
            new_node["sequence"] += n2i["sequence"]
            new_node["name"] = "merged_%s-%s" % (n1i["name"], n2i["name"])

        all_edges = list(edges_by_name.values())
        new_name = new_node["name"]

        for v in all_edges:
            if v["to"] == n1:
                v["to"] = new_name
            if v["from"] == n1:
                v["from"] = new_name
            if v["to"] == n2:
                v["to"] = new_name
            if v["from"] == n2:
                v["from"] = new_name

        all_edges = [e for e in all_edges if e["from"] != e["to"]]

        del nodes_by_name[n1]
        del nodes_by_name[n2]
        nodes_by_name[new_name] = new_node

        fadj = defaultdict(list)
        radj = defaultdict(list)
        for e in all_edges:
            fadj[e["from"]].append(e["to"])
            radj[e["to"]].append(e["from"])
        edges_by_name = {e["from"] + "_" + e["to"]: e for e in all_edges}
        nodes_to_check = list([n for n in nodes_by_name.values()
                               if len(fadj[n["name"]]) == 1 and len(radj[fadj[n["name"]][0]]) == 1])

    all_nodes = sorted(nodes_by_name.values(), key=lambda x: x["topo_order"])
    all_edges = sorted(edges_by_name.values(), key=lambda y: (y["from"], y["to"]))

    for n in all_nodes:
        del n["topo_order"]

    return list(all_nodes), list(all_edges)


def cleanup_and_add_paths(all_nodes, all_edges):
    """
    Cleanup and final check of graph, add paths
    :param all_nodes: list of nodes
    :param all_edges: list of edges
    :return: nodes, edges, path_infos
    """

    all_nodes, all_edges = contract_paths(all_nodes, all_edges)

    all_edges = sorted(all_edges, key=lambda e_: (e_["from"], e_["to"]))

    fadj = defaultdict(list)
    radj = defaultdict(list)

    for e in all_edges:
        fadj[e["from"]].append(e["to"])
        radj[e["to"]].append(e["from"])

    sequences_to_check = set()
    for n in all_nodes:
        if "mark" in n:
            del n["mark"]
        if "sequences" in n:
            n["sequences"] = list(n["sequences"])
            if not n["sequences"]:
                del n["sequences"]
            else:
                for s in n["sequences"]:
                    sequences_to_check.add(s)

    for e in all_edges:
        if "sequences" in e:
            e["sequences"] = list(e["sequences"])
            for s in e["sequences"]:
                sequences_to_check.add(s)
            if not e["sequences"]:
                del e["sequences"]

    nodes_by_name = {n["name"]: n for n in all_nodes}
    edges_by_name = {e["from"] + "_" + e["to"]: e for e in all_edges}

    def has_seq_edge(node_, sequence_):
        for n_ in fadj[node_]:
            en = node_ + "_" + n_
            if en in edges_by_name and "sequences" in edges_by_name[en] and sequence_ in edges_by_name[en]["sequences"]:
                return True
        for n_ in radj[node_]:
            en = n_ + "_" + node_
            if en in edges_by_name and "sequences" in edges_by_name[en] and sequence_ in edges_by_name[en]["sequences"]:
                return True
        return False

    path_infos = []
    nodes_covered_by_sequences = set()
    for s in sequences_to_check:
        nodes_with_this_sequence = [n for n in all_nodes
                                    if "sequences" in n and s in n["sequences"] or has_seq_edge(n["name"], s)]

        def can_go(n1, n2, sq):
            return ("sequences" not in nodes_by_name[n1] or sq in nodes_by_name[n1]["sequences"] or list(nodes_by_name[n1]["sequences"]) == ["REF"]) and \
                   ("sequences" not in nodes_by_name[n2] or sq in nodes_by_name[n2]["sequences"] or list(nodes_by_name[n2]["sequences"]) == ["REF"]) and \
                   ("sequences" not in edges_by_name[n1 + "_" + n2] or
                    sq in edges_by_name[n1 + "_" + n2]["sequences"]) and \
                   ("sequences" in nodes_by_name[n1] or
                    "sequences" in nodes_by_name[n2] or
                    "sequences" in edges_by_name[n1 + "_" + n2])

        def node_sequence_length(n_):
            if "reference" in nodes_by_name[n_]:
                _, start_, end_ = parse_region(nodes_by_name[n_]["reference"])
                return end_ - start_ + 1
            else:
                return len(nodes_by_name[n_]["sequence"])

        path_count = 1
        while nodes_with_this_sequence:
            start = nodes_with_this_sequence.pop(0)

            current_name = start["name"]
            while radj[current_name]:
                preds = [r for r in radj[current_name] if can_go(r, current_name, s)]
                if not preds:
                    break
                elif len(preds) > 1:
                    raise Exception("Sequence %s is ambiguous before node %s" % (s, current_name))
                else:
                    current_name = preds[0]

            path = [current_name]
            while fadj[current_name]:
                succs = [r for r in fadj[current_name] if can_go(current_name, r, s)]
                if not succs:
                    break
                elif len(succs) > 1:
                    raise Exception("Sequence %s is ambiguous after node %s" % (s, current_name))
                else:
                    current_name = succs[0]
                    path.append(current_name)

            for node in path:
                nodes_covered_by_sequences.add(node)

            nodes_with_this_sequence = [n for n in nodes_with_this_sequence
                                        if n["name"] not in set(path)]
            path_infos.append({
                "sequence": s,
                "nodes": path,
                "nucleotide_length": sum([node_sequence_length(n) for n in path
                                          if n not in {"source", "sink"}]),
                "path_id": s + "|" + str(path_count)
            })
            path_count += 1

    for x in path_infos:
        if x["nodes"][0] in fadj["source"]:
            x["nodes"].insert(0, "source")
            nodes_covered_by_sequences.add("source")
        if x["nodes"][-1] in radj["sink"]:
            x["nodes"].append("sink")
            nodes_covered_by_sequences.add("sink")

    path_infos = sorted(path_infos, key=lambda pi: pi["path_id"])

    for n in all_nodes:
        if "sequences" in n:
            n["sequences"] = sorted(n["sequences"])

    for e in all_edges:
        if "sequences" in e:
            e["sequences"] = sorted(e["sequences"])

    all_node_names = set(nodes_by_name.keys())
    if all_node_names - nodes_covered_by_sequences:
        raise Exception("Some nodes that aren't covered by any sequences: %s" %
                        str(all_node_names - nodes_covered_by_sequences))

    return all_nodes, all_edges, path_infos
