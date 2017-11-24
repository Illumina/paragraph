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
# Topological sort for DAG
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

from collections import defaultdict


def topological_sort(all_nodes,
                     all_edges,
                     source_name="source",
                     sink_name="sink"):
    """
    Topological sort of graph (necessary for use with GSSW), also adds source and sink if necessary
    :param all_nodes: list of nodes
    :param all_edges: list of edges
    :param source_name: name of source node
    :param sink_name: name of sink node
    :return: sorted_nodes, sorted_edges -- first element in sorted_nodes is source, last is sink
    """

    if not all_nodes:
        return all_nodes, all_edges

    all_edges = sorted(all_edges, key=lambda e_: (e_["from"], e_["to"]))

    fadj = defaultdict(list)
    radj = defaultdict(list)
    for e in all_edges:
        fadj[e["from"]].append(e["to"])
        radj[e["to"]].append(e["from"])

    ancestor_set = {source_name, sink_name}

    def fix_ancestor(anode, nodes_seen=None):
        if not nodes_seen:
            nodes_seen = []
        if anode in nodes_seen and anode not in ancestor_set:
            raise Exception("Graph has a cycle at %s (%s)" % (anode, str(nodes_seen)))
        nodes_seen.append(anode)
        if anode in ancestor_set:
            return
        else:
            if not radj[anode]:
                all_edges.insert(0, {"from": source_name, "to": anode})
                radj[anode].append(source_name)
                fadj[source_name].append(anode)
            else:
                ancestor_set.add(anode)
                for prev in radj[anode]:
                    if prev not in ancestor_set:
                        fix_ancestor(prev, nodes_seen)

    for n in all_nodes:
        fix_ancestor(n["name"])

    terminal_set = {source_name, sink_name}

    def fix_terminal(tnode, nodes_seen=None):
        if not nodes_seen:
            nodes_seen = []
        if tnode in nodes_seen and tnode not in ancestor_set:
            raise Exception("Graph has a cycle at %s (%s)" % (tnode, str(nodes_seen)))
        nodes_seen.append(tnode)
        if tnode in terminal_set:
            return
        else:
            if not fadj[tnode]:
                all_edges.insert(0, {"from": tnode, "to": sink_name})
                fadj[tnode].append(sink_name)
                radj[sink_name].append(tnode)
            else:
                terminal_set.add(tnode)
                for succ in fadj[tnode]:
                    if succ not in terminal_set:
                        fix_terminal(succ, nodes_seen)

    for n in all_nodes:
        fix_terminal(n["name"])

    nodes_by_name = {n["name"]: n for n in all_nodes}
    if source_name not in nodes_by_name:
        all_nodes.insert(0, {
            "name": source_name,
            "sequence": "N" * 10
        })
        nodes_by_name[source_name] = all_nodes[0]

    if sink_name not in nodes_by_name:
        all_nodes.append({
            "name": sink_name,
            "sequence": "N" * 10
        })
        nodes_by_name[sink_name] = all_nodes[-1]

    ordered_nodes = []

    def visit(start):
        nodes_by_name[start]["mark"] = "t"
        for x in sorted(fadj[start]):
            if "mark" not in nodes_by_name[x]:
                visit(x)
            elif "mark" in nodes_by_name[x] and nodes_by_name[x]["mark"] == "t":
                raise Exception("Graph has a cycle at %s -> %s" % (start, x))
        nodes_by_name[start]["mark"] = "p"
        ordered_nodes.insert(0, start)

    visit(all_nodes[0]["name"])

    all_nodes = [nodes_by_name[nn] for nn in ordered_nodes]
    for n in all_nodes:
        del n["mark"]

    return all_nodes, all_edges
