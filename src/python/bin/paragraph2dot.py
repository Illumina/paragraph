#!/usr/bin/env python3
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
#
# June 2017
#
# Convert paragraph JSON file to graphviz / dot file
#
# Usage:
#
# For usage instructions run with option --help
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import json
import argparse


def main():
    parser = argparse.ArgumentParser("paragraph2dot.py")

    parser.add_argument("input", help="Input paragraph JSON file", nargs=1)
    parser.add_argument("output", help="Output dot file", nargs=1)

    args = parser.parse_args()

    with open(args.input[0]) as f:
        paragraph_data = json.load(f)

    output_file = open(args.output[0], "w")

    print("digraph paragraph_export {", file=output_file)

    node_mapping = {}
    for i, node in enumerate(paragraph_data["nodes"]):
        if "reference" in node:
            node_shape = "ellipse"
            node_label = str(i) + ": " + node["reference"]
        else:
            node_shape = "box"
            node_sequence = node["sequence"]
            if len(node_sequence) > 30:
                node_sequence = node_sequence[:14] + "..." + node_sequence[-14:]
            node_label = str(i) + ": " + node_sequence

        if "sequences" in node:
            node_label += " (%s)" % str(node["sequences"])

        node_mapping[node["name"]] = "node_%i" % i

        print("node_%i [label=\"%s\" shape=%s];" % (i, node_label, node_shape),
              file=output_file)

    for e in paragraph_data["edges"]:
        edgelabel = ""
        if "sequences" in e:
            edgelabel += str(e["sequences"])
        print("%s -> %s  [label=\"%s\"];" % (node_mapping[e["from"]],
                                             node_mapping[e["to"]],
                                             edgelabel),
              file=output_file)

    print("}", file=output_file)


if __name__ == "__main__":
    main()
