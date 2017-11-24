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
# October 2017
#
# Run grmpy genotyper on multiple paragraph output
# Output will be one JSON for one sample with genotyping information
#
# Use multiparagraph.py as template
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
import gzip
import itertools
import json
import logging
import multiprocessing
import os
import pipes
import subprocess
import tempfile
import traceback
import copy

import findgrm  # pylint: disable=unused-import
from grm.helpers import LoggingWriter
from grm.graph_typing import sample_manifest


def run_single_sample_grmpy(event_and_args):
    """
    Run grmpy for one single sample on all variants using multiple threads
    """
    event = event_and_args[0]
    args = event_and_args[1]
    sample_name = args.sample_name
    tempfiles = []
    exception = ""
    error_log = ""
    error = False
    gt_result = {}
    try:
        # construct grmpy input
        current_json = {"ID": event[args.id_identifier],
                        "nodes": event["graph"]["nodes"],
                        "edges": event["graph"]["edges"],
                        "target_regions": event["graph"]["target_regions"],
                        "sequencenames": event["graph"]["sequencenames"],
                        "samples": {},
                        "eventinfo": {}}
        current_json["samples"][sample_name] = {}
        if "vcf_records" in event["graph"]:
            current_json["eventinfo"]["vcf_records"] = event["graph"]["vcf_records"]
        for e in ["chrom", "start", "end", "type", "n_ev"]:
            try:
                current_json["eventinfo"][e] = event[e]
            except KeyError:
                current_json["eventinfo"][e] = None

        try:
            current_json["samples"][sample_name]["read_counts_by_edge"] = event["graph"]["read_counts_by_edge"]
        except KeyError:
            current_json["samples"][sample_name]["read_counts_by_edge"] = {}

        tf = tempfile.NamedTemporaryFile(dir=args.scratch_dir, mode="wt", suffix=".json", delete=False)
        tempfiles.append(tf.name)
        json.dump(current_json, tf)
        tf.close()

        # grmpy output
        output_name = tf.name + ".output.json"
        tempfiles.append(output_name)
        error_log = tf.name + ".output.log"
        tempfiles.append(error_log)

        commandline = args.grmpy
        commandline += " -r %s" % pipes.quote(args.reference)
        commandline += " -m %s" % pipes.quote(args.manifest)
        commandline += " -p %s" % pipes.quote(tf.name)
        commandline += " -o %s" % pipes.quote(output_name)
        commandline += " --genotype-error-rate %f" % args.genotype_error_rate
        commandline += " --min-overlap-bases %i" % args.min_overlap_bases
        commandline += " --max-read-times %i" % args.max_read_times
        if args.use_em:
            commandline += " --useEM"
        commandline += " --log-file %s" % pipes.quote(error_log)

        o = subprocess.check_output(commandline, shell=True, stderr=subprocess.STDOUT)

        try:
            o = o.decode("utf-8")
        except:  # pylint: disable=bare-except
            o = str(o)

        for line in o.split("\n"):
            if line:
                logging.warning(line)

        with open(output_name, "rt") as f:
            gt_result = json.load(f)
            f.close()

    except Exception:  # pylint: disable=broad-except
        logging.error("Exception when running grmpy on %s", str(event))
        logging.error('-' * 60)
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        logging.error('-' * 60)
        exception = traceback.format_exc()
        error = True
    except BaseException:  # pylint: disable=broad-except
        logging.error("Exception when running grmpy on %s", str(event))
        logging.error('-' * 60)
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        logging.error('-' * 60)
        exception = traceback.format_exc()
        error = True
    finally:
        if error:
            try:
                with open(error_log, "rt") as f:
                    for line in f:
                        line = str(line).strip()
                        logging.error(line)
            except:  # pylint: disable=bare-except
                pass
            gt_result["error"] = {
                "exception": exception
            }
        if not args.keep_scratch:
            for f in tempfiles:
                try:
                    os.remove(f)
                except:  # pylint: disable=bare-except
                    pass
    return gt_result


def make_argument_parser():
    """
    :return: an argument parser
    """
    parser = argparse.ArgumentParser("Multigrmpy.py")

    parser.add_argument("-m", "--manifest", help="paragraph output json manifest. Format: sample_name, json_path.",
                        type=str, dest="manifest", required=True)

    parser.add_argument("-o", "--output", help="Output directory name",
                        type=str, dest="output", required=True)

    parser.add_argument("-r", "--reference-sequence", help="Reference FASTA",
                        type=str, dest="reference", required=True)

    parser.add_argument("--id_identifier", help="ID key identifier in Json.",
                        type=str, dest="id_identifier", default="ID")

    parser.add_argument("--event-threads", "-t", dest="threads", type=int, default=multiprocessing.cpu_count(),
                        help="Number of events to process in parallel.")

    parser.add_argument("--keep-scratch", dest="keep_scratch", default=None, action="store_true",
                        help="Do not delete temp files.")

    parser.add_argument("--scratch-dir", dest="scratch_dir", default=None,
                        help="Directory for temp files")

    parser.add_argument("--grmpy", dest="grmpy", default=os.path.join(findgrm.GRM_BASE, "bin", "grmpy"),
                        type=str, help="Path to the grmpy executable")

    parser.add_argument("--logfile", dest="logfile", default=None,
                        help="Write logging information into file rather than to stderr")

    verbosity_options = parser.add_mutually_exclusive_group(required=False)

    verbosity_options.add_argument("--verbose", dest="verbose", default=False, action="store_true",
                                   help="Raise logging level from warning to info.")

    verbosity_options.add_argument("--quiet", dest="quiet", default=False, action="store_true",
                                   help="Set logging level to output errors only.")

    stat_options = parser.add_mutually_exclusive_group(required=False)

    stat_options.add_argument("--genotype-error-rate", dest="genotype_error_rate", default=0.01,
                              type=float, help="Fixed genotype error rate for breakpoint genotyping.")

    stat_options.add_argument("--min-overlap-bases", dest="min_overlap_bases", default=16, type=int,
                              help="Minimum overlap bases used in estimating poisson model parameters.")

    stat_options.add_argument("--max-read-times", dest="max_read_times", default=40, type=int,
                              help="Max times of total reads in one sample for a breakpoint. Multiplied by depth.")
    stat_options.add_argument("--useEM", dest="use_em", default=False,
                              action="store_true", help="Use EM for genotyping.")

    return parser


def run(raw_args):
    """ Run wrapper """
    if raw_args.verbose:
        loglevel = logging.INFO
    elif raw_args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.WARNING

    # reinitialize logging
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=raw_args.logfile, format='%(asctime)s %(levelname)-8s %(message)s', level=loglevel)

    try:
        manifest = sample_manifest.load_manifest(raw_args.manifest)
    except:  # pylint: disable=bare-except
        logging.error("Error in loading paragraph manifest.")
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        raise
    if not manifest:
        raise Exception("Empty paragraph manifest.")

    os.makedirs(raw_args.output, exist_ok=True)
    for element in manifest:
        args = copy.copy(raw_args)

        tmp_manifest = tempfile.NamedTemporaryFile(
            dir=args.scratch_dir, mode="wt", suffix=".manifest", delete=False)
        tmp_manifest.write(element.to_string() + "\n")
        tmp_manifest.close()
        args.manifest = tmp_manifest.name

        paragraph_json_path = element.path

        logging.info("Loading paragraph output for sample " + element.name)
        if paragraph_json_path.endswith(".gz"):
            paragraph_json_file = gzip.open(paragraph_json_path, "rt")
        else:
            paragraph_json_file = open(paragraph_json_path, "rt")
        paragraph_variants = json.load(paragraph_json_file)
        logging.info("Loaded paragraph output Json.")
        paragraph_json_file.close()

        args.sample_name = element.name
        args.output = os.path.join(
            raw_args.output, element.name + ".genotype.json.gz")
        logging.info("Number of events: %i", len(paragraph_variants))

        with multiprocessing.Pool(args.threads) as pool:
            results = pool.map(run_single_sample_grmpy, zip(
                paragraph_variants, itertools.repeat(args)))

        logging.info("Finished genotyping. Merging output...")
        if args.output.endswith(".gz"):
            of = gzip.open(args.output, "wt")
        else:
            of = open(args.output, "wt")

        json.dump(results, of, sort_keys=True,
                  indent=4, separators=(',', ': '))
        of.close()
        if not args.keep_scratch:
            try:
                os.remove(tmp_manifest.name)
            except:  # pylint: disable=bare-except
                pass
        logging.info("Done genotyping of " + element.name)


def main():
    parser = make_argument_parser()
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
