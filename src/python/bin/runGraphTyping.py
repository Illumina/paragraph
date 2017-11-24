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
# September 2017
#
# Main function to run paragraph and grmpy from given JSON/VCF and bam file list
#
# Basic usage:
#  runGraphTyping.py -i <input json or vcf> -b <bam manifest> -r <reference file> -o <working directory>
#
# For usage instructions run with option --help
#
# Author:
#
# Sai Chen <schen6@illumina.com>
#

import logging
import argparse
import os
import sys
import json
import traceback
import importlib
import findgrm  # pylint: disable=unused-import

from grm.helpers import LoggingWriter
from grm.vcf2paragraph import convert_vcf_to_json
from grm.graph_typing import idx_depth_wrapper
from grm.graph_typing import sample_manifest


class LoggingContext(object):  # pylint: disable=too-few-public-methods
    """ Helper context to handle our other bin submodules messing with our logging config """

    def __init__(self, logging_target=None):
        self.logging_target = logging_target

    def __enter__(self):
        logging.shutdown()
        importlib.reload(logging)

    def __exit__(self, *args):
        logging.shutdown()
        importlib.reload(logging)
        logging.basicConfig(filename=self.logging_target,
                            format='%(asctime)s: %(message)s',
                            level=logging.INFO)


def main():
    parameters = load_parameters()

    bam_manifest = sample_manifest.load_manifest(parameters.manifest, parameters.relative_manifest_path)

    if not bam_manifest:
        raise Exception("Empty BAM list. Please check your option -b")

    # convert vcf to Json if necessary
    extension = os.path.splitext(parameters.input)[1]
    if extension == ".vcf" or extension == ".vcf.gz":
        logging.info("Input is a vcf. Column3 will be recognized as ID.")
        json_input_path = os.path.join(parameters.output, "variants.json")
        if os.path.isfile(json_input_path):
            logging.info("Json variant representations generated at:\n" + json_input_path + "\nUse this json instead.")
        else:
            try:
                event_list = convert_vcf_to_json(parameters.input,
                                                 parameters.reference,
                                                 read_length=parameters.read_length,
                                                 max_ref_node_length=parameters.max_ref_node_length,
                                                 graph_type=parameters.graph_type,
                                                 split_type=parameters.split_type,
                                                 retrieve_ref_sequence=parameters.retrieve_reference_sequence,
                                                 threads=parameters.threads)
                with open(json_input_path, "wt") as json_file:
                    json.dump(event_list, json_file,
                              sort_keys=True,
                              indent=4,
                              separators=(',', ': '))
            except Exception:  # pylint: disable=W0703
                logging.error("VCF conversion failed.")
                traceback.print_exc(file=LoggingWriter(logging.ERROR))
                raise
            logging.info("Done. Graph Json stored at: \n" + json_input_path)
        parameters.input = json_input_path
    elif extension != ".json" and extension != ".json.gz":
        raise Exception("Unknown input file extension %s for %s" %
                        (extension, parameters.input))

    # prepare idxdepth for genotyping if this information is not included in manifest
    if sample_manifest.count_missing_stats(bam_manifest) == len(bam_manifest):
        logging.info("Missing depth and read length info in manifest.")
        new_manifest_path = os.path.join(parameters.output, "idx_updated.manifest")
        if os.path.isfile(new_manifest_path):
            logging.info("Depth and read length info generated at:\n" +
                         new_manifest_path + "\nUse this manifest instead.")
        else:
            logging.info("Generating idxdepth info for each sample. This may take a while...")
            idx_depth_dir = os.path.join(parameters.output, "idxdepth")
            os.makedirs(idx_depth_dir, exist_ok=True)
            try:
                idx_depth_wrapper.generate_idxdepth_json(
                    bam_manifest,
                    parameters.reference,
                    idx_depth_dir,
                    threads=parameters.threads,
                    idxdepth_binary=parameters.idxdepth_binary)
            except Exception:  # pylint: disable=W0703
                logging.error("Idxdepth failed for at least one sample.")
                traceback.print_exc(file=LoggingWriter(logging.ERROR))
                raise
            logging.info("Done generating idxdepth info.")
            sample_manifest.generate_and_update_manifest(idx_depth_dir, bam_manifest, new_manifest_path)
            logging.info("Done. New manifest stored at: \n" + new_manifest_path)
        parameters.manifest = new_manifest_path
        bam_manifest = sample_manifest.load_manifest(parameters.manifest, parameters.relative_manifest_path)
    if sample_manifest.count_missing_stats(bam_manifest) > 0:
        raise Exception(
            "Missing either depth or read length in one sample. Please check: " + parameters.manifest)

    # run multi paragraph to for read counts
    logging.info("Running multi-paragraph for graph read counts...")
    import multiparagraph
    for bam_element in bam_manifest:
        logging.info("Processing sample %s...", bam_element.name)
        os.makedirs(os.path.join(
            parameters.output, "paragraph", bam_element.name))
        try:
            parser = multiparagraph.make_argument_parser()
            argv = list(map(str, [parameters.input,
                                  "-b", bam_element.path,
                                  "-r", parameters.reference,
                                  "-o", os.path.join(parameters.output, "paragraph",
                                                     bam_element.name, "raw_pg.json"),
                                  "-t", parameters.threads,
                                  "--paragraph", parameters.paragraph_binary,
                                  "--logfile", os.path.join(parameters.output, "paragraph", bam_element.name,
                                                            "multiparagraph.log")]))
            if parameters.output_everything:
                argv.append("-E")
            if parameters.paragraph_binary:
                argv.extend(["--paragraph", parameters.paragraph_binary])
            if parameters.max_events:
                argv.extend(["--max-events", parameters.max_events])
            if parameters.min_length:
                argv.extend(["--min-length", parameters.min_length])
            args = parser.parse_args(argv)
            with LoggingContext(parameters.logfile):
                multiparagraph.run(args)
        except Exception as e:  # pylint: disable=W0703
            logging.error(
                "Multi paragraph failed for sample " + bam_element.name)
            traceback.print_exc(file=LoggingWriter(logging.ERROR))
            raise e
    logging.info("Done running multi-paragraph.")

    # generate manifest
    pg_manifest_path = os.path.join(
        parameters.output, "paragraph/paragraph.manifest")
    with open(pg_manifest_path, "wt") as pg_manifest_file:
        for bam_element in bam_manifest:
            pg_path = os.path.join(
                parameters.output, "paragraph", bam_element.name, "raw_pg.json")
            pg_element = bam_element
            pg_element.path = pg_path
            pg_manifest_file.write(bam_element.to_string() + "\n")

    # run multigrmpy for genotyping
    logging.info("Running multi-grmpy for genotyping...")
    import multigrmpy
    genotyping_logfile = os.path.join(parameters.output, "multigrmpy.log")
    try:
        parser = multigrmpy.make_argument_parser()
        argv = list(map(str, ["-m", pg_manifest_path,
                              "-r", parameters.reference,
                              "-o", parameters.output,
                              "-t", parameters.threads,
                              "--id_identifier", parameters.id_identifier,
                              "--genotype-error-rate", parameters.genotype_error_rate,
                              "--min-overlap-bases", parameters.min_overlap_bases,
                              "--max-read-times", parameters.max_read_times,
                              "--logfile", parameters.logfile]))
        if parameters.grmpy_binary:
            argv.extend(["--grmpy", parameters.grmpy_binary])
        if parameters.use_em:
            argv.append("--useEM")
        args = parser.parse_args(argv)
        with LoggingContext(parameters.logfile):
            multigrmpy.run(args)
    except Exception as e:  # pylint: disable=W0703
        logging.error("Genotyping failed. See %s", genotyping_logfile)
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        raise e
    logging.info(
        "Done running genotyping. Please check output at: " + parameters.output)


def load_parameters():
    prog_name = sys.argv[0]
    parser = argparse.ArgumentParser(description=prog_name)
    parser.add_argument("-i", "--input", help="Input file of graphs. Can be either JSON or VCF.",
                        type=str, dest="input", required=True)
    parser.add_argument("-o", "--out", help="Output directory. Will create one if it doesn't exist.",
                        type=str, dest="output", required=True)
    parser.add_argument("-m", "--manifest", help="List of bam files. Should have 2 coluns with name and path.",
                        type=str, dest="manifest", required=True)
    parser.add_argument("-r", "--reference-sequence",
                        help="Reference FASTA for checking REF and resolving <DEL>",
                        type=str, dest="reference", required=True)

    paragraphs = parser.add_argument_group("Paragraph options")
    paragraphs.add_argument("-E", "--output-everything",
                            dest="output_everything", default=False, action="store_true",
                            help="Output full alignment information. WARNING: this produces very large output files.")
    paragraphs.add_argument("--max-events", dest="max_events", type=int,
                            default=None, help="Only do the first n events.")
    paragraphs.add_argument("--min-length", dest="min_length", type=int, default=0, help="Minimum event length.")

    grmpys = parser.add_argument_group("Grmpy options")
    grmpys.add_argument("--id_identifier", help="ID key identifier in Json.",
                        type=str, dest="id_identifier", default="ID")
    grmpys.add_argument("--genotype-error-rate", dest="genotype_error_rate", default=0.01,
                        type=float, help="Fixed genotype error rate for breakpoint genotyping.")
    grmpys.add_argument("--min-overlap-bases", dest="min_overlap_bases", default=16, type=int,
                        help="Minimum overlap bases used in estimating poisson model parameters.")
    grmpys.add_argument("--max-read-times", dest="max_read_times", default=40, type=int,
                        help="Max times of total reads in one sample for a breakpoint. Multiplied by depth.")
    grmpys.add_argument("--useEM", dest="use_em", default=False, action="store_true", help="Use EM for genotyping.")

    v2p = parser.add_argument_group("Options for making graph from vcf")
    v2p.add_argument("--graph-type", choices=["alleles", "haplotypes"], default="alleles", dest="graph_type",
                     help="Select the type of complex graph to generate. Same as the option --graph-type in vcf2paragraph.")
    v2p.add_argument("--vcf-split", default="lines", dest="split_type", choices=["lines", "full", "by_id"],
                     help="Mode for splitting the input VCF: lines (default) -- one graph per record ;"
                     " full -- one graph for the whole VCF ; by_id -- use the VCF id column to group adjacent records")
    v2p.add_argument("--retrieve-reference-sequence", help="Retrieve reference sequence for REF nodes",
                     action="store_true", dest="retrieve_reference_sequence", default=False)
    v2p.add_argument("-l", "--max-ref-node-length", dest="max_ref_node_length", type=int, default=1000,
                     help="Maximum length of reference nodes before they get padded and truncated.")
    v2p.add_argument("-p", "--read-length", dest="read_length", type=int, default=150,
                     help="Read length -- this can be used to add reference padding for disambiguation.")

    binaries = parser.add_argument_group("Custom binaries")
    binaries.add_argument("--grmpy-binary", dest="grmpy_binary", default=os.path.join(findgrm.GRM_BASE, "bin", "grmpy"),
                          type=str, help="Path to the grmpy executable")
    binaries.add_argument("--paragraph-binary", dest="paragraph_binary",
                          default=os.path.join(findgrm.GRM_BASE, "bin", "paragraph"), help="Path to the paragraph executable")
    binaries.add_argument("--idxdepth-binary", dest="idxdepth_binary", default=os.path.join(findgrm.GRM_BASE, "bin", "idxdepth"),
                          help="Path to the idxdepth executable")

    miscs = parser.add_argument_group("Misc")
    miscs.add_argument("--threads", dest="threads", default=1, type=int, help="Number of threads to use.")
    miscs.add_argument("--relative-manifest-path", dest="relative_manifest_path", action="store_true",
                       default=False, help="Add manifest directory to paths in manifest.")
    miscs.add_argument("--logfile", dest="logfile", default=None, help="Specify log file")
    miscs.add_argument("--keep-scratch", dest="keep_scratch", default=None,
                       action="store_true", help="Do not delete temp files.")

    args = parser.parse_args()

    parameter_encoding = ' '.join(sys.argv[1:])

    os.makedirs(args.output, exist_ok=True)

    if not args.logfile:
        args.logfile = os.path.join(args.output, "GraphTyping.log")

    logging.shutdown()
    importlib.reload(logging)
    logging.basicConfig(filename=args.logfile, format='%(asctime)s: %(message)s', level=logging.INFO)

    logging.info("Starting with these parameters: %s", parameter_encoding)
    return args


if __name__ == "__main__":
    logging.basicConfig(format='%(asctime)s: %(message)s', level=logging.INFO)
    try:
        main()
    except Exception as e:  # pylint: disable=W0703
        print(str(e), file=sys.stderr)
        logging.error(str(e))
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        exit(1)
