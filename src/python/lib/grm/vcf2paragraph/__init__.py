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
# June 2017
#
# Wrapper for VCF to graph conversion
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com> and Sai Chen <schen6@illumina.com>
#

import os
import multiprocessing
import logging
import itertools
import traceback
import hashlib
import tempfile

import pysam.bcftools  # pylint: disable=E0401
import pysam  # pylint: disable=E0401

from grm.vcfgraph import VCFGraph, NoVCFRecordsException
from grm.vcfgraph import graphUtils
from grm.vcfgraph.graphContainer import GraphContainer
from grm.helpers import parse_region
from grm.helpers import LoggingWriter


def add_reference_information(paragraph_dict, reference_fasta):
    """ Adds reference sequence information to reference nodes """
    fasta = pysam.FastaFile(reference_fasta)
    for n in paragraph_dict["nodes"]:
        if "reference" in n:
            chrom, start, end = parse_region(n["reference"])
            n["reference_sequence"] = fasta.fetch(chrom, start - 1, end).upper()


def convert_vcf(vcf,
                ref,
                ins_info_key,
                target_regions=None,
                ref_node_padding=150,
                ref_node_max_length=1000,
                allele_graph=False,
                simplify=True,
                alt_paths=False,
                alt_splitting=False):
    """
    Convert a single VCF file to a graph dictionary
    :param vcf: file name of the VCF file
    :param ref: reference FASTA file name
    :param target_regions: target region list
    :param ref_node_padding: padding / read length
    :param ref_node_max_length: maximum length before splitting a reference node
    :param allele_graph: add edges between any compatible allele pair, not just haplotypes from input
    :param simplify: simplify the graph
    :param alt_paths: Add all possible non-reference paths to the graph
    :param alt_splitting: also split long alt nodes (e.g. long insertions)
    :return: dictionary containing JSON graph
    """
    graph = GraphContainer("Graph from %s" % vcf)
    indexed_vcf = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf.gz")
    try:
        indexed_vcf.close()
        # noinspection PyUnresolvedReferences
        pysam.bcftools.view(vcf, "-o", indexed_vcf.name, "-O", "z", catch_stdout=False)  # pylint: disable=no-member
        # noinspection PyUnresolvedReferences
        pysam.bcftools.index(indexed_vcf.name)  # pylint: disable=no-member

        regions = map(parse_region, target_regions) if target_regions else [(None,)*3]
        for (chrom, start, end) in regions:
            if chrom is not None:
                logging.info("Starting work on region: %s:%d-%d", chrom, start, end)
            try:
                vcfGraph = VCFGraph.create_from_vcf(
                    ref, indexed_vcf.name, ins_info_key, chrom, start, end, ref_node_padding, allele_graph)
            except NoVCFRecordsException:
                logging.info("Region %s:%d-%d has no VCF records, skipping.", chrom, start, end)
                continue
            logging.info("Constructed VCF graph:\n%s", str(vcfGraph))
            chromGraph = vcfGraph.get_graph(allele_graph)
            if ref_node_max_length:
                graphUtils.split_ref_nodes(chromGraph, ref_node_max_length, ref_node_padding)
                if alt_splitting:
                    graphUtils.split_alt_nodes(chromGraph, ref_node_max_length, ref_node_padding)

            if simplify:
                graphUtils.remove_empty_nodes(chromGraph)
                graphUtils.combine_nodes(chromGraph)
                # Disable edge label simplification for now. May use node-label short-cut later
                # graphUtils.remove_redundant_edge_labels(graph)
            chromGraph.check()

            graphUtils.add_graph(graph, chromGraph)
    finally:
        os.remove(indexed_vcf.name)

    graph.target_regions = target_regions or graph.get_reference_regions()
    graphUtils.add_source_sink(graph)
    graphUtils.add_ref_path(graph)
    if alt_paths:
        graphUtils.add_alt_paths(graph)
    graph.check()
    return graph.json_dict()


def convert_vcf_to_json(args, alt_paths=False):
    """ Convert VCF file to list of JSON graphs
    This function converts each record separately (as opposed to all records jointly like convert_vcf)
    :param vcf_name: VCF file name
    :param reference: reference FASTA name
    :param read_length: read length parameter for VCF conversion
    :param max_ref_node_length: length before splitting ref nodes
    :param graph_type: type of graph for vcf2paragraph
    :param split_type: "lines", "full", "by_id", "superloci"
    :param retrieve_ref_sequence: retrieve the reference sequence
    :param alt_splitting: split long alt nodes
    :param threads: number of processes to use
    :param alt_paths: generate paths for ALT alleles
    :return header, vcf_records, events -- vcf header, matched lists of VCF records and an event
    """

    header, records, block_ids = parse_vcf_lines(args.input, args.read_length, args.split_type)
    variants = [None] * len(block_ids)

    to_process = []
    try:
        for record_block in records:
            tf = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf.gz")
            tf.close()
            to_process.append(tf.name)
            with pysam.VariantFile(tf.name, 'w', header=header) as vcf_out:
                for record in record_block:
                    vcf_out.write(record)

        params = {"reference": args.reference,
                  "read_length": args.read_length,
                  "max_ref_node_length": args.max_ref_node_length,
                  "graph_type": args.graph_type,
                  "retrieve_reference_sequence": args.retrieve_reference_sequence,
                  "alt_splitting": args.alt_splitting,
                  "alt_paths": alt_paths,
                  "ins_info_key": args.ins_info_key}

        with multiprocessing.Pool(args.threads) as pool:
            variants = pool.map(run_vcf2paragraph, zip(to_process, itertools.repeat(params)))

        if any([x is None for x in variants]):
            raise Exception("Event conversion failed for at least one VCF record.")

        for v, b in zip(variants, block_ids):
            v["ID"] = b
    finally:
        for p in to_process:
            try:
                os.remove(p)
            except:  # pylint: disable=bare-except
                pass

    return header, records, variants


def parse_vcf_lines(vcf_path, read_length=150, split_type="full"):
    """
    Split VCF into parts

    Store header and records separately
    :param vcf_path: path to VCF file
    :param read_length: length of reads = minimum padding length
    :param split_type: "lines", "full", "by_id", "superloci"

    :return: vcf header, block records, block IDs
    """
    sha = hashlib.sha256()
    blocksize = 65536
    with open(vcf_path, 'rb') as vcf_file:
        file_buffer = vcf_file.read(blocksize)
        while file_buffer:
            sha.update(file_buffer)
            file_buffer = vcf_file.read(blocksize)
    vcf_id = os.path.basename(vcf_path) + "@" + str(sha.hexdigest())

    vcf_file = pysam.VariantFile(vcf_path)
    vcf_header = vcf_file.header
    # pylint: disable=no-member
    if 'GRMPY_ID' not in vcf_header.info:
        vcf_header.add_line('##INFO=<ID=GRMPY_ID,Number=1,Type=String,Description="Graph ID '
                            'for linking to genotypes.json.gz; matches record.graphinfo.ID in there.">')
    # pylint: enable=no-member

    records = []
    block_ids = []
    prev_id = ""
    current_chr = None
    previous_end = None

    for record in vcf_file:
        if record.pos < read_length:
            raise Exception("Distance between vcf position and chrom start is smaller than read length.")
        if split_type == "full":
            bid = vcf_id + ":0"
            record.info["GRMPY_ID"] = bid
            if not records:
                records = [[record]]
                block_ids.append(bid)
            else:
                records[0].append(record)
        elif split_type == "lines":
            bid = vcf_id + ":" + str(len(records) + 1)
            record.info["GRMPY_ID"] = bid
            records.append([record])
            block_ids.append(bid)
        elif split_type == "by_id":
            if not record.id:
                # records without ID are just appended separately
                # make sure we add these separately
                bid = vcf_id + ":" + str(len(records) + 1)
                record.info["GRMPY_ID"] = bid
                records.append([record])
                block_ids.append(bid)
                prev_id = None
            elif record.id == prev_id:
                bid = block_ids[-1]
                record.info["GRMPY_ID"] = bid
                records[-1].append(record)
            else:
                bid = vcf_id + ":" + str(len(records) + 1)
                record.info["GRMPY_ID"] = bid
                records.append([record])
                block_ids.append(bid)
                prev_id = record.id
        elif split_type == "superloci":
            if not current_chr or record.chrom != current_chr or \
               not previous_end or record.pos > previous_end + read_length:
                bid = vcf_id + ":" + str(len(records) + 1)
                record.info["GRMPY_ID"] = bid
                records.append([record])
                block_ids.append(bid)
                logging.debug("New superlocus started at %s:%i -- previous : %s:%d ; read_length = %i",
                              str(record.chrom), record.pos, str(current_chr), previous_end, read_length)
            else:
                bid = block_ids[-1]
                record.info["GRMPY_ID"] = bid
                records[-1].append(record)
            current_chr = record.chrom
            previous_end = record.stop
            if not previous_end or previous_end < record.pos:
                previous_end = record.pos
        else:
            raise Exception("Unknown VCF splitting type: %s" % split_type)
    vcf_file.close()
    return vcf_header, records, block_ids


def run_vcf2paragraph(event_and_args):
    """
    run vcf2paragraph for one single variant
    """
    event = event_and_args[0]
    params = event_and_args[1]
    tempfiles = []
    result = {}

    try:
        logging.debug("Converting: %s", str(event))

        result["graph"] = convert_vcf(
            event,
            params["reference"],
            params["ins_info_key"],
            None,
            ref_node_padding=params["read_length"],
            ref_node_max_length=params["max_ref_node_length"],
            allele_graph=params["graph_type"] == "alleles",
            alt_splitting=params["alt_splitting"],
            alt_paths=params["alt_paths"])

        chrom = None
        start = None
        end = None
        if "vcf_records" in result["graph"]:
            for r in result["graph"]["vcf_records"]:
                if chrom is None:
                    chrom = r["chrom"]
                else:
                    assert chrom == r["chrom"]
                if start is None:
                    start = r["pos"]
                else:
                    start = min(start, r["pos"])
                if end is None:
                    end = r["end"]
                else:
                    end = max(end, r["end"])
        else:
            for tr in result["graph"]["target_regions"]:
                c, s, e = parse_region(tr)
                if chrom is None:
                    chrom = c
                else:
                    assert chrom == c
                if start is None:
                    start = s
                else:
                    start = min(start, s)
                if end is None:
                    end = e
                else:
                    end = max(end, e)

        assert chrom is not None
        assert start is not None
        assert end is not None
        result["chrom"] = chrom
        result["start"] = start
        result["end"] = end
    except Exception:  # pylint: disable=broad-except
        logging.error("Exception when running vcf2paragraph on %s", str(event))
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        raise
    finally:
        for x in tempfiles:
            try:
                os.remove(x)
            except:  # pylint: disable=bare-except
                pass

    if params["retrieve_reference_sequence"]:
        add_reference_information(result["graph"], params["reference"])

    return result
