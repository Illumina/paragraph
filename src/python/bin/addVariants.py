#!/usr/bin/env python3

# coding=utf-8
#
# Copyright (c) 2016-2019 Illumina, Inc.
# All rights reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# You may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied
# See the License for the specific language governing permissions and limitations
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
