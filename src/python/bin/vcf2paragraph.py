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
# Convert VCF to paragraph JSON file
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

# noinspection PyUnresolvedReferences
import findgrm  # pylint: disable=unused-import

from grm.vcf2paragraph import convert_haploid_vcf, convert_allele_vcf, add_reference_information


def make_argument_parser():
    """
    :return: one argument parser
    """

    parser = argparse.ArgumentParser("vcf2paragraph.py")

    parser.add_argument("input", help="Input VCF / BCF file", nargs=1)
    parser.add_argument("output", help="Output JSON file", nargs=1)

    parser.add_argument("-r", "--reference-sequence", help="Reference FASTA for checking REF and resolving <DEL>",
                        type=str, dest="ref", required=True)

    common = parser.add_argument_group("Common VCF graph options")
    common.add_argument("-g", "--graph-type", choices=["alleles", "haplotypes"], default="haplotypes",
                        dest="graph_type", help="Select the type of graph to generate")
    common.add_argument("-R", "--retrieve-reference-sequence", help="Retrieve reference sequence for REF nodes",
                        action="store_true", dest="retrieve_reference_sequence", default=False)
    common.add_argument("-l", "--max-ref-node-length", dest="max_ref_len", type=int, default=1000,
                        help="Maximum length of reference nodes before they get padded and truncated.")
    common.add_argument("-p", "--read-length", dest="read_len", type=int, default=150,
                        help="Read length -- this can be used to add reference padding for disambiguation.")
    common.add_argument("-T", "--target-region", dest="target_regions", default=[], action="append",
                        help="Target regions for read retrieval; also, reference nodes will not be clipped "
                             "inside target regions.")

    haploid = parser.add_argument_group("Haploid VCF graph options")
    haploid.add_argument("--ref-path", help="Add edges for the reference path.",
                         action="store_true", dest="ref_paths", default=False)
    haploid.add_argument("--crossovers", help="Add crossover edges (connect directly adjacent"
                                              "nodes even when no sample gives these combinations).",
                         action="store_true", dest="crossovers", default=False)

    return parser


def run(args):
    if args.graph_type == "alleles":
        output = convert_allele_vcf(args.input[0],
                                    args.target_regions,
                                    args.read_len,
                                    args.max_ref_len)
    else:
        output = convert_haploid_vcf(args.input[0],
                                     args.ref,
                                     args.target_regions,
                                     args.read_len,
                                     args.max_ref_len,
                                     args.ref_paths,
                                     args.crossovers)
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
