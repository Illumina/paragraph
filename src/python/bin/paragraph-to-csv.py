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
# January 2018
#
# Extract variant genotypes from paragraph JSON and output as csv format
#
# Usage:
#   python3 paragraph-to-csv.py <paragraph json> [options]
#
# For usage instructions run with option --help
#
# Author:
#
# Sai Chen <schen6@illumina.com>
#

# noinspection PyUnresolvedReferences

import argparse
import json
import gzip
import tempfile
import os


def main():
    parser = make_argument_parser()
    args = parser.parse_args()
    run(args)


def make_argument_parser():
    parser = argparse.ArgumentParser("paragraph-to-vcf.py")
    parser.add_argument("input", nargs="+")
    parser.add_argument("-o", "--output", type=str, dest="output",
                        help="Output csv file name. If not specified, will be printed to standard output.", )
    parser.add_argument("--genotype-only", dest="genotype_only", default=False,
                        action="store_true", help="Only output genotypes for each variant.")
    return parser


def run(args):
    """
    Run converter
    """
    if args.input[0].endswith(".gz"):
        paragraph_json_file = gzip.open(args.input[0], 'rt')
    else:
        paragraph_json_file = open(args.input[0], 'rt')
    paragraph_json = json.load(paragraph_json_file)
    paragraph_json_file.close()

    if not paragraph_json:
        raise Exception("Fail to load paragraph JSON.")

    if args.genotype_only:
        annotation = "#FORMAT=GT"
    else:
        annotation = "#FORMAT=GT:FILTER:DP"

    samples = []
    header = "#ID"
    for event in paragraph_json:
        if "samples" not in event:
            continue
        if "population" in event:
            header += "\tHWE\tCallRate\tPassRate\tMinorAF"
        for sample in event["samples"]:
            header += "\t" + sample
            samples.append(sample)
        break
    if not samples:
        raise Exception("No valid results in paragraph JSON")

    output_lines = []
    for event in paragraph_json:
        if "graphinfo" not in event:
            output_lines.append("ERROR")
            continue
        line = ""
        if "ID" in event["graphinfo"]:
            line = event["graphinfo"]["ID"]
        if not line:
            line = "."

        if "population" in event:
            pop_json = event["population"]
            pop_stat = []
            for k in ["hwe", "call_rate"]:
                if k in pop_json:
                    pop_stat.append(pop_json[k])
                else:
                    pop_stat.append("NA")
            # calculate pass rate
            num_pass = 0
            num_samples = 0
            for sample in samples:
                num_samples += 1
                sample_stat = event["samples"][sample]["gt"]["GT"]
                if "filter" in sample_stat:
                    if sample_stat["filter"] == "PASS":
                        num_pass += 1
            if num_samples:
                pop_stat.append(num_pass / num_samples)
            else:
                pop_stat.append("NA")
            # AF
            if "allele_frequencies" in pop_json:
                if pop_json["allele_frequencies"]:
                    pop_stat.append(min(pop_json["allele_frequencies"]))
            else:
                pop_stat.append("NA")
            line += "\t" + '\t'.join(map(str, pop_stat))

        for sample in samples:
            gt_json = event["samples"][sample]["gt"]
            sample_stat = [gt_json["GT"]]
            if not args.genotype_only:
                for k in ["filter", "num_reads"]:
                    if k in gt_json:
                        sample_stat.append(gt_json[k])
                    else:
                        sample_stat.append(".")
            line += "\t" + ':'.join(map(str, sample_stat))
        output_lines.append(line)

    if args.output:
        output_file = open(args.output, 'wt')
    else:
        output_file = tempfile.NamedTemporaryFile(mode="wt", suffix=".csv", delete=False)

    output_file.write(annotation + "\n")
    output_file.write(header + "\n")
    for line in output_lines:
        output_file.write(line + "\n")
    output_file.close()

    if not args.output:
        os.system("cat " + output_file.name)
        os.system("rm -f " + output_file.name)


if __name__ == '__main__':
    main()
