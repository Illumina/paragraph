#!/usr/bin/env python3
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
#
# April 2018
#
# Usage:
#
# For usage instructions run with option --help
#
# Author:
#
# Felix Schlesinger <fschlesinger@illumina.com>
#

import json
import argparse
import logging

# noinspection PyUnresolvedReferences
import findgrm  # pylint: disable=unused-import
from grm.vcfgraph import graphUtils, variants
from grm.helpers import load_json


def make_argument_parser():
    parser = argparse.ArgumentParser(description="Add new variants into a graph as nodes/edges")
    parser.add_argument("graph", help="Input graph JSON file")
    parser.add_argument("output", help="Output JSON file", type=argparse.FileType('w'))
    parser.add_argument("--variants", default=None,
                        help="JSON file with variant calls to use instead of calls in input graph")
    parser.add_argument("-v", "--verbose", action="count", default=0,
                        help="More logging; May be given twice for even more logging")
    return parser


def run(args):
    levels = [logging.WARNING, logging.INFO, logging.DEBUG]
    logging.basicConfig(format='%(message)s', level=levels[min(args.verbose, len(levels)-1)])

    graphDict = load_json(args.graph)
    graph = graphUtils.load_json(graphDict)
    if args.variants:
        varJson = load_json(args.variants)
        if "variants" not in varJson:
            raise Exception("No variants in variant JSON")
        varDict = varJson["variants"]
    else:
        varDict = graphDict.get("variants", {})
        if not varDict:
            logging.warning("No variants in graph")
            print(varDict.keys())
    variants.add_variants(graph, varDict)
    graphUtils.remove_empty_nodes(graph)
    json.dump(graph.json_dict(), args.output, sort_keys=True)


def main():
    args = make_argument_parser().parse_args()
    run(args)


if __name__ == "__main__":
    main()
