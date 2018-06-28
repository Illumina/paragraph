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
# Run paragraph tools on multiple sites on multiple samples
#   This script is suitable for a small list of sites and samples.
#   To run ParaGRAPH on large-scale sites and samples, please refer to the Snakemake template at doc/multi-samples.md
# Output will be one JSON for all samples with genotyping information
#
# Usage:
#
# For usage instructions run with option --help
#
# Author:
#
# Sai Chen <schen6@illumina.com>
#

import argparse
import logging
import traceback

import findgrm  # pylint: disable=unused-import
from grm.helpers import LoggingWriter
from grm.vcfgraph import vcfupdate


def make_argument_parser():
    """
    :return: an argument parser
    """
    parser = argparse.ArgumentParser("Multigrmpy.py")

    parser.add_argument("-i", "--input", help="Input VCF file of variants.",
                        type=str, dest="input", required=True)

    parser.add_argument("-g", "--grmpy", help="JSON output from multigrmpy.py / grmpy",
                        type=str, dest="input_grm", required=True)

    parser.add_argument("-o", "--output", help="Output file name.", type=str, dest="output", default="-")

    parser.add_argument("--logfile", dest="logfile", default=None,
                        help="Write logging information into file rather than to stderr")

    verbosity_options = parser.add_mutually_exclusive_group(required=False)

    verbosity_options.add_argument("--verbose", dest="verbose", default=False, action="store_true",
                                   help="Raise logging level from warning to info.")

    verbosity_options.add_argument("--quiet", dest="quiet", default=False, action="store_true",
                                   help="Set logging level to output errors only.")

    verbosity_options.add_argument("--debug", dest="debug", default=False, action="store_true",
                                   help="Log debug level events.")

    return parser


def run(args):
    """
    :run the wrapper
    """
    if args.verbose:
        loglevel = logging.INFO
    elif args.quiet:
        loglevel = logging.ERROR
    elif args.debug:
        loglevel = logging.DEBUG
    else:
        loglevel = logging.WARNING

    # reinitialize logging
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=args.logfile, format='%(asctime)s %(levelname)-8s %(message)s', level=loglevel)

    try:
        grmpyOutput = vcfupdate.read_grmpy(args.input_grm)
        vcfupdate.update_vcf_from_grmpy(args.input, grmpyOutput, args.output)
    except Exception:  # pylint: disable=W0703
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        raise


def main():
    parser = make_argument_parser()
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
