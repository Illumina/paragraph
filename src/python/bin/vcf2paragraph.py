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

from grm.vcf2paragraph import convert_vcf, add_reference_information


def make_argument_parser():
    """
    :return: one argument parser
    """

    parser = argparse.ArgumentParser("vcf2paragraph.py")

    parser.add_argument("input", help="Input VCF / BCF file", nargs=1)
    parser.add_argument("output", help="Output JSON file", nargs=1)

    parser.add_argument("-r", "--reference-sequence", type=str, dest="ref", required=True,
                        help="Reference FASTA for checking REF and resolving <DEL>")
    parser.add_argument("-v", "--verbose", action="count", default=0,
                        help="More logging; May be given twice for even more logging.")
    common = parser.add_argument_group("Common VCF graph options")
    common.add_argument("-g", "--graph-type", choices=["alleles", "haplotypes"],
                        default="haplotypes", dest="graph_type",
                        help="Select the type of graph to generate.")
    common.add_argument("-R", "--retrieve-reference-sequence", action="store_true",
                        dest="retrieve_reference_sequence", default=False,
                        help="Retrieve reference sequence for REF nodes")
    common.add_argument("-l", "--max-ref-node-length", dest="max_ref_len", type=int, default=1000,
                        help="Maximum length of reference nodes before they get padded and truncated.")
    common.add_argument("-p", "--read-length", dest="read_len", type=int, default=150,
                        help="Read length -- this can be used to add reference padding for disambiguation.")
    common.add_argument("-T", "--target-region", dest="target_regions", default=[], action="append",
                        help="Target regions for read retrieval")
    common.add_argument("--ins-info-key", dest="ins_info_key", default="SEQ",
                        type=str, help="Key for symbolic <INS> in INFO field")
    common.add_argument("--alt-paths", dest="alt_paths", default=False, action="store_true",
                        help="Create all possible ALT paths in addition to reference paths.")
    common.add_argument("--alt-splitting", dest="alt_splitting", default=False, action="store_true",
                        help="Also split long alternate alleles (e.g. long insertions)")
    common.add_argument("--recursion-limit", dest="recursion_limit", default=None, type=int,
                        help="Set the recursion limit ( O(expected number of nodes of the graph) for large graphs"
                             " -- this is required for sorting )")
    return parser


def run(args):
    levels = [logging.WARNING, logging.INFO, logging.DEBUG]
    logging.basicConfig(format='%(message)s', level=levels[min(args.verbose, len(levels) - 1)])

    if args.recursion_limit:
        import sys
        sys.setrecursionlimit(args.recursion_limit)

    output = convert_vcf(args.input[0],
                         args.ref,
                         args.ins_info_key,
                         args.target_regions,
                         args.read_len,
                         args.max_ref_len,
                         args.graph_type == "alleles",
                         alt_paths=args.alt_paths,
                         alt_splitting=args.alt_splitting)
    if args.retrieve_reference_sequence:
        add_reference_information(output, args.ref)

    with open(args.output[0], "w") as output_file:
        json.dump(output, output_file, sort_keys=True)


def main():
    parser = make_argument_parser()
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
