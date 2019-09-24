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
import argparse
import gzip
import json
import logging
import multiprocessing
import os
import pipes
import subprocess
import tempfile
import traceback
import re
import pysam

import findgrm  # pylint: disable=unused-import
from grm.helpers import LoggingWriter
from grm.vcf2paragraph import convert_vcf_to_json
from grm.graph_templates import make_graph
from grm.vcfgraph import vcfupdate


def load_graph_description(args):
    """
    load graph description from either vcf or json
    """
    event_list = []
    extension = os.path.splitext(args.input)[1]
    if extension == ".gz":
        file_type = os.path.splitext(os.path.splitext(args.input)[0])[1]
        extension = file_type + ".gz"
    if extension == ".vcf" or extension == ".vcf.gz":
        logging.info("Input is a vcf. Converting to JSON with graph description...")
        try:
            converted_json_path = os.path.join(args.output, "variants.json.gz")
            header, records, event_list = convert_vcf_to_json(args, alt_paths=True)

            vcf_with_event_ids_path = os.path.join(args.output, "variants.vcf.gz")
            logging.info("Saving: %s.", vcf_with_event_ids_path)
            with pysam.VariantFile(vcf_with_event_ids_path, 'w', header=header) as vcf_with_event_ids_file:
                for record in [record for record_block in records for record in record_block]:
                    vcf_with_event_ids_file.write(record)

            logging.info("Saving: %s.", converted_json_path)
            with gzip.open(converted_json_path, "wt") as converted_json_file:
                json.dump(event_list, converted_json_file, sort_keys=True, indent=4, separators=(',', ': '))
        except Exception:  # pylint: disable=W0703
            logging.error("VCF to JSON conversion failed.")
            traceback.print_exc(file=LoggingWriter(logging.ERROR))
            raise
        logging.info("Done. Graph Json stored at %s", converted_json_path)
    elif extension == ".json" or extension == ".json.gz":
        if extension == ".json":
            json_file = open(args.input, 'r')
        else:
            json_file = gzip.open(args.input, 'r')
        event_list = json.load(json_file)
        num_converted_event = 0
        # if JSON has no graph description
        for event in event_list:
            if "graph" not in event:
                if "nodes" not in event and "edges" not in event:
                    try:
                        event["type"], event["graph"] = make_graph(args.reference, event)
                    except Exception:  # pylint: disable=W0703
                        logging.error("Failed to make graph for JSON event.")
                        traceback.print_exc(file=LoggingWriter(logging.ERROR))
                        raise
                    num_converted_event += 1
        if num_converted_event:
            logging.info("Constructed graph for %d events in JSON.", num_converted_event)
        json_file.close()
    else:
        raise Exception("Unknown input file extension %s for %s. Only VCF or JSON is allowed!" %
                        (extension, args.input))

    tempfiles = []
    logging.info("Saving %d graph json files", len(event_list))
    graph_id = 0
    for event in event_list:
        input_json_file = tempfile.NamedTemporaryFile(dir=args.scratch_dir, mode="wt", suffix=".json", delete=False)
        tempfiles.append(input_json_file.name)
        if "graph" in event:
            graph = event["graph"]
            if "ID" not in graph or not graph["ID"]:
                if "ID" in event:
                    graph["ID"] = event["ID"]
                else:
                    graph["ID"] = os.path.basename(args.input) + ":" + str(graph_id)

            # graph id gives index into original events s.t. we can find them even if
            # only some of them have an ID already
            graph_id += 1
            json.dump(graph, input_json_file, indent=4, separators=(',', ': '))
        else:
            json.dump(event, input_json_file, indent=4, separators=(',', ': '))
        input_json_file.close()
    return tempfiles


def make_argument_parser():
    """
    :return: an argument parser
    """
    parser = argparse.ArgumentParser("Multigrmpy.py")

    parser.add_argument("-i", "--input", help="Input file of variants. Must be either JSON or VCF.",
                        type=str, dest="input", required=True)

    parser.add_argument("-m", "--manifest", help="Manifest of samples with path and bam stats.",
                        type=str, dest="manifest", required=True)

    parser.add_argument("-o", "--output", help="Output directory.", type=str, dest="output", required=True)

    parser.add_argument("-A", "--write-alignments", dest="write_alignments",
                        help="Write alignment JSON files into the output folder (large!).",
                        default=False, action="store_true")

    parser.add_argument("--infer-read-haplotypes", dest="infer_read_haplotypes",
                        help="Infer read haplotype paths",
                        default=False, action="store_true")

    parser.add_argument("-r", "--reference-sequence", help="Reference genome fasta file.",
                        type=str, dest="reference", required=True)

    parser.add_argument("--threads", "-t", dest="threads", type=int, default=multiprocessing.cpu_count(),
                        help="Number of events to process in parallel.")

    parser.add_argument("--keep-scratch", dest="keep_scratch", default=None, action="store_true",
                        help="Do not delete temp files.")

    parser.add_argument("--scratch-dir", dest="scratch_dir", default=None,
                        help="Directory for temp files")

    parser.add_argument("--grmpy", dest="grmpy", default=os.path.join(findgrm.GRM_BASE, "bin", "grmpy"),
                        type=str, help="Path to the grmpy executable")

    parser.add_argument("--logfile", dest="logfile", default=None,
                        help="Write logging information into file rather than to stderr")

    parser.add_argument("--graph-sequence-matching", dest="graph_sequence_matching",
                        default=False, help="Use graph aligner.")

    parser.add_argument("--klib-sequence-matching", dest="klib_sequence_matching",
                        default=False, help="Use klib smith waterman aligner.")

    parser.add_argument("--kmer-sequence-matching", dest="kmer_sequence_matching",
                        default=False, help="Use kmer aligner.")

    parser.add_argument("--bad-align-uniq-kmer-len", dest="bad_align_uniq_kmer_len", default=0,
                        help="Kmer length for uniqueness check during read filtering.")

    parser.add_argument("--no-alt-splitting", dest="alt_splitting", default=True, action="store_false",
                        help="Keep long insertion sequences in the graph rather than trimming them at the read / padding length.")

    verbosity_options = parser.add_mutually_exclusive_group(required=False)

    verbosity_options.add_argument("--verbose", dest="verbose", default=False, action="store_true",
                                   help="Raise logging level from warning to info.")

    verbosity_options.add_argument("--quiet", dest="quiet", default=False, action="store_true",
                                   help="Set logging level to output errors only.")

    verbosity_options.add_argument("--debug", dest="debug", default=False, action="store_true",
                                   help="Log debug level events.")

    stat_options = parser.add_argument_group('gt-parameter-group')

    stat_options.add_argument("-G", "--genotyping-parameters", dest="genotyping_parameters", default="",
                              type=str, help="JSON string or file with genotyping model parameters.")

    stat_options.add_argument("-M", "--max-reads-per-event", dest="max_reads_per_event", default=0,
                              type=int, help="Maximum number of reads to process for a single event.")

    vcf2json_options = parser.add_argument_group('vcf2json-option-group')

    vcf2json_options.add_argument("--vcf-split", default="lines", dest="split_type", choices=["lines", "full", "by_id", "superloci"],
                                  help="Mode for splitting the input VCF: lines (default) -- one graph per record ;"
                                  " full -- one graph for the whole VCF ;"
                                  " by_id -- use the VCF id column to group adjacent records ;"
                                  " superloci -- split VCF into blocks where records are separated by at least read-length characters")
    vcf2json_options.add_argument("-p", "--read-length", dest="read_length", default=150, type=int,
                                  help="Read length -- this is used to add reference padding when converting VCF to graphs.")

    vcf2json_options.add_argument("-l", "--max-ref-node-length", dest="max_ref_node_length", type=int, default=300,
                                  help="Maximum length of reference nodes before they get padded and truncated.")

    vcf2json_options.add_argument("--retrieve-reference-sequence", help="Retrieve reference sequence for REF nodes",
                                  action="store_true", dest="retrieve_reference_sequence", default=False)

    vcf2json_options.add_argument("--graph-type", choices=["alleles", "haplotypes"], default="alleles", dest="graph_type",
                                  help="Type of complex graph to generate. Same as --graph-type in vcf2paragraph.")
    vcf2json_options.add_argument("--ins-info-key", type=str, default="SEQ", dest="ins_info_key",
                                  help="Key in INFO field to indicate sequence of symbolic <INS>")

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

    os.makedirs(args.output, exist_ok=True)
    if args.scratch_dir:
        os.makedirs(args.scratch_dir, exist_ok=True)

    # reinitialize logging
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=args.logfile, format='%(asctime)s %(levelname)-8s %(message)s', level=loglevel)

    # manifest sanity check
    with open(args.manifest) as manifest_file:
        headers = {"id": False, "path": False, "idxdepth": False, "depth": False,
                   "read length": False, "sex": False, "depth variance": False, "depth sd": False}
        for line in manifest_file:
            if line.startswith("#"):
                line = line[1:]
            line = line.strip()
            fields = re.split('\t|,', line)
            for field in fields:
                if field not in headers:
                    header_str = ",".join(headers)
                    raise Exception("Illegal header name %s. Allowed headers:\n%s" % (field, header_str))
                headers[field] = True
            if not headers["id"] or not headers["path"]:
                raise Exception("Missing header \"id\" or \"path\" in manifest")
            if not headers["idxdepth"]:
                if not headers["depth"] or not headers["read length"]:
                    raise Exception("Missing header \"idxdepth\", or \"depth\" and \"read length\" in manifest.")
            break

    # prepare input graph description
    result_json_path = os.path.join(args.output, "genotypes.json.gz")
    try:
        graph_files = load_graph_description(args)
        commandline = " -r %s" % pipes.quote(args.reference)
        commandline += " -m %s" % pipes.quote(args.manifest)
        commandline += " -o %s" % pipes.quote(result_json_path)
        commandline += " -z"

        if args.genotyping_parameters:
            commandline += " -G %s" % pipes.quote(args.genotyping_parameters)
        if args.max_reads_per_event:
            commandline += " -M %s" % pipes.quote(str(args.max_reads_per_event))
        if args.threads >= 1:
            commandline += " -t %s" % pipes.quote(str(args.threads))
        if args.graph_sequence_matching:
            commandline += " --graph-sequence-matching %s" % pipes.quote(str(args.graph_sequence_matching))
        if args.klib_sequence_matching:
            commandline += " --klib-sequence-matching %s" % pipes.quote(str(args.klib_sequence_matching))
        if args.kmer_sequence_matching:
            commandline += " --kmer-sequence-matching %s" % pipes.quote(str(args.kmer_sequence_matching))
        if int(args.bad_align_uniq_kmer_len):
            commandline += " --bad-align-uniq-kmer-len %s" % pipes.quote(str(args.bad_align_uniq_kmer_len))
        if args.write_alignments:
            alignment_directory = os.path.join(args.output, "alignments")
            os.makedirs(alignment_directory, exist_ok=True)
            if not os.path.isdir(alignment_directory):
                raise Exception("Cannot create alignment output directory: {}".format(alignment_directory))
            commandline += " --alignment-output-folder !%s" % pipes.quote(alignment_directory)
        if args.infer_read_haplotypes:
            commandline += " --infer-read-haplotypes"

        if args.verbose:
            commandline += " --log-level=info"
        elif args.quiet:
            commandline += " --log-level=error"
        elif args.debug:
            commandline += " --log-level=debug"
        else:
            commandline += " --log-level=warning"

        grmpy_log = pipes.quote(os.path.join(args.output, "grmpy.log"))
        commandline += " --log-file " + grmpy_log
        commandline += " --log-async no"

        commandline += " -g"
        for graph in graph_files:
            commandline += "\n%s" % pipes.quote(graph)

        response_file = tempfile.NamedTemporaryFile(dir=args.scratch_dir, mode="wt", suffix=".txt", delete=False)
        response_file.write(commandline)
        response_file.flush()

        commandline = args.grmpy + " --response-file=%s" % pipes.quote(response_file.name)

        logging.info("Starting: %s", commandline)

        subprocess.check_call(commandline, shell=True, stderr=subprocess.STDOUT)

    except Exception:  # pylint: disable=W0703
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        raise

    try:
        sample_names = []
        with open(args.manifest) as manifest_file:
            id_index = -1
            for line in manifest_file:
                line = line.rstrip()
                if line.startswith('#'):
                    line = line[1:]
                fields = re.split('\t|,', line)
                if id_index == -1:
                    id_index = fields.index("id")
                    continue
                sample_names.append(fields[id_index])
        if args.input.endswith("vcf") or args.input.endswith("vcf.gz"):
            grmpyOutput = vcfupdate.read_grmpy(result_json_path)
            result_vcf_path = os.path.join(args.output, "genotypes.vcf.gz")
            vcf_input_path = os.path.join(args.output, "variants.vcf.gz")
            if not os.path.exists(vcf_input_path) or not os.path.isfile(vcf_input_path):
                vcf_input_path = args.input
            vcfupdate.update_vcf_from_grmpy(vcf_input_path, grmpyOutput, result_vcf_path, sample_names)
    except Exception:  # pylint: disable=W0703
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        raise


def main():
    parser = make_argument_parser()
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
