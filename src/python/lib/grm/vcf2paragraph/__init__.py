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
from collections import defaultdict
import tempfile
import multiprocessing
import json
import logging
import itertools
import gzip
import traceback

import intervaltree
import pysam.bcftools
import pysam

from grm.vcfgraph import HaploidVCFGraph, AlleleGraph
from grm.vcfgraph.topological_sort import topological_sort
from grm.vcfgraph.ref_node_split import split_long_reference_nodes
from grm.vcfgraph.cleanup import cleanup_and_add_paths
from grm.helpers import parse_region
from grm.helpers import LoggingWriter


def add_reference_information(paragraph_dict, reference_fasta):
    """ Adds reference sequence information to reference nodes """
    fasta = pysam.FastaFile(reference_fasta)
    for n in paragraph_dict["nodes"]:
        if "reference" in n:
            chrom, start, end = parse_region(n["reference"])
            n["reference_sequence"] = fasta.fetch(chrom, start - 1, end).upper()


def convert_allele_vcf(vcf_path,
                       target_regions=None,
                       ref_node_padding=150,
                       ref_node_max_length=1000):
    """
    Convert a single VCF file to a graph dictionary
    :param vcf_path: file name of the VCF file
    :param target_regions: target region list
    :param ref_node_padding: padding / read length
    :param ref_node_max_length: maximum length before splitting a reference node
    :return: dictionary containing JSON graph
    """
    if not target_regions:
        target_regions = []
    output = {"model_name": "Allele graph from %s" % vcf_path,
              "sequencenames": [],
              "nodes": [],
              "edges": [],
              "vcf_records": []}

    indexed_vcf = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf.gz")
    try:
        indexed_vcf.close()
        # noinspection PyUnresolvedReferences
        pysam.bcftools.view(vcf_path, "-o", indexed_vcf.name, "-O", "z", catch_stdout=False)  # pylint: disable=no-member
        # noinspection PyUnresolvedReferences
        pysam.bcftools.index(indexed_vcf.name)  # pylint: disable=no-member

        all_nodes = []
        all_edges = []
        regions_added = set()
        region_list = list(map(parse_region, target_regions))
        if not list(region_list):
            region_list = [(None, None, None)]

        for i, (chrom, start, end) in enumerate(region_list):
            graph = AlleleGraph.create_from_vcf(
                indexed_vcf.name, chrom, start, end, ref_node_padding)
            nodes, edges = graph.get_nodes_and_edges()
            output["sequencenames"] += list(graph.all_sequences)
            output["vcf_records"] += graph.vcflist
            if ref_node_max_length:
                nodes, edges = split_long_reference_nodes(
                    nodes, edges, ref_node_max_length, ref_node_padding)
            nodes, edges = topological_sort(nodes, edges)
            all_nodes += nodes
            all_edges += edges
            if nodes:
                regions_added.add(i)
        if len(regions_added) > 1:
            output["nodes"], output["edges"] = topological_sort(all_nodes, all_edges)
        else:
            output["nodes"], output["edges"] = all_nodes, all_edges
        output["sequencenames"] = list(set(output["sequencenames"]))
    finally:
        os.remove(indexed_vcf.name)

    # determine target regions from all reference segments
    if not target_regions:
        reference_trees = defaultdict(intervaltree.IntervalTree)
        for n in all_nodes:
            if "reference" in n:
                n_chr, n_start, n_end = parse_region(n["reference"])
                reference_trees[n_chr].addi(n_start, n_end + 1)
        target_regions = []
        for chrom, tree in reference_trees.items():
            tree.merge_overlaps()
            for iv in tree.items():
                target_regions.append("%s:%i-%i" %
                                      (chrom, iv.begin, iv.end - 1))

    output["nodes"], output["edges"], output["paths"] = cleanup_and_add_paths(output["nodes"], output["edges"])
    output["target_regions"] = sorted(target_regions)
    output["sequencenames"] = sorted(output["sequencenames"])

    return output


def convert_haploid_vcf(vcf_path,
                        ref,
                        target_regions=None,
                        ref_node_padding=150,
                        ref_node_max_length=1000,
                        ref_paths=False,
                        crossovers=False):
    """
    Convert a single VCF file to a graph dictionary
    :param vcf_path: file name of the VCF file
    :param ref: reference FASTA file name
    :param target_regions: target region list
    :param ref_node_padding: padding / read length
    :param ref_node_max_length: maximum length before splitting a reference node
    :param ref_paths: add reference paths even if not in VCF
    :param crossovers: add crossover paths
    :return: dictionary containing JSON graph
    """
    if not target_regions:
        target_regions = []
    output = {"model_name": "Graph from %s" % vcf_path,
              "sequencenames": [],
              "nodes": [],
              "edges": []}

    indexed_vcf = tempfile.NamedTemporaryFile(delete=False, suffix=".vcf.gz")
    try:
        indexed_vcf.close()
        # noinspection PyUnresolvedReferences
        pysam.bcftools.view(vcf_path, "-o", indexed_vcf.name, "-O", "z", catch_stdout=False)  # pylint: disable=no-member
        # noinspection PyUnresolvedReferences
        pysam.bcftools.index(indexed_vcf.name)  # pylint: disable=no-member

        all_nodes = []
        all_edges = []
        regions_added = set()
        region_list = list(map(parse_region, target_regions))
        if not list(region_list):
            region_list = [(None, None, None)]
        for i, (chrom, start, end) in enumerate(region_list):
            graph = HaploidVCFGraph.create_from_vcf(
                ref, indexed_vcf.name, chrom, start, end, ref_node_padding)
            nodes, edges = graph.get_nodes_and_edges(ref_paths, crossovers)
            output["sequencenames"] += list(graph.all_sequences)
            if ref_node_max_length:
                nodes, edges = split_long_reference_nodes(
                    nodes, edges, ref_node_max_length, ref_node_padding)
            nodes, edges = topological_sort(nodes, edges)
            all_nodes += nodes
            all_edges += edges
            if nodes:
                regions_added.add(i)
        if len(regions_added) > 1:
            output["nodes"], output["edges"] = topological_sort(
                all_nodes, all_edges)
        else:
            output["nodes"], output["edges"] = all_nodes, all_edges
        output["sequencenames"] = list(set(output["sequencenames"]))
    finally:
        os.remove(indexed_vcf.name)

    # determine target regions from all reference segments
    if not target_regions:
        reference_trees = defaultdict(intervaltree.IntervalTree)
        for n in all_nodes:
            if "reference" in n:
                n_chr, n_start, n_end = parse_region(n["reference"])
                reference_trees[n_chr].addi(n_start, n_end + 1)
        target_regions = []
        for chrom, tree in reference_trees.items():
            tree.merge_overlaps()
            for iv in tree.items():
                target_regions.append("%s:%i-%i" %
                                      (chrom, iv.begin, iv.end - 1))

    output["nodes"], output["edges"], output["paths"] = cleanup_and_add_paths(
        output["nodes"], output["edges"])
    output["target_regions"] = sorted(target_regions)
    output["sequencenames"] = sorted(output["sequencenames"])

    return output


def convert_vcf_to_json(vcf_name,
                        reference,
                        read_length=150,
                        max_ref_node_length=300,
                        graph_type="alleles",
                        split_type="lines",
                        retrieve_ref_sequence=False,
                        threads=1):
    """ Convert VCF file to list of JSON graphs
    This function converts each record separately (as opposed to all records jointly like convert_vcf)
    :param vcf_name: VCF file name
    :param reference: reference FASTA name
    :param read_length: read length parameter for VCF conversion
    :param max_ref_node_length: length before splitting ref nodes
    :param graph_type: type of graph for vcf2paragraph
    :param split_type: "lines", "full", "by_id"
    :param retrieve_ref_sequence: retrieve the reference sequence
    :param threads: number of processes to use
    """

    header, records = parse_vcf_lines(vcf_name, read_length, split_type)

    params = {"header": header,
              "reference": reference,
              "read_length": read_length,
              "max_ref_node_length": max_ref_node_length,
              "graph_type": graph_type,
              "retrieve_reference_sequence": retrieve_ref_sequence}

    with multiprocessing.Pool(threads) as pool:
        variants = pool.map(run_vcf2paragraph, zip(records, itertools.repeat(params)))

    if any([x is None for x in variants]):
        raise Exception("Event conversion failed for at least one VCF record.")

    return variants


def parse_vcf_lines(vcf_path, read_length=150, split_type="full"):
    """
    Split VCF into parts

    Store header and records separately
    :param vcf_path: path to VCF file
    :param split_type: "lines", "full", "by_id"
    """
    vcf_file = pysam.VariantFile(vcf_path)
    header = str(vcf_file.header)
    records = []
    prev_id = ""
    for record in vcf_file.fetch(None, None, None):
        if record.pos < read_length:
            raise Exception("Distance between vcf position and chrom start is smaller than read length.")
        if split_type != "lines":
            if split_type == "full":
                if not records:
                    records[0] = str(record)
                else:
                    records[0] += str(record)
            else:
                if not record.id:
                    # records without ID are just appended separately
                    # make sure we add these separately
                    records.append(str(record))
                    prev_id = None
                elif record.id == prev_id:
                    records[-1] += "\n" + str(record)
                else:
                    records.append(str(record))
                    prev_id = record.id
        elif split_type == "lines":
            records.append(str(record))
        else:
            raise Exception("Unknown VCF splitting type: %s" % split_type)
    vcf_file.close()
    return header, records


def run_vcf2paragraph(event_and_args):
    """
    run vcf2paragraph for one single variant
    """
    event = event_and_args[0]
    params = event_and_args[1]
    tempfiles = []
    result = {}

    try:
        logging.info(event)

        # prepare input
        tf = tempfile.NamedTemporaryFile(
            mode="wt", suffix=".json", delete=False)
        tempfiles.append(tf.name)
        tf.write(params["header"])
        tf.write(event)
        tf.close()

        # output
        if params["graph_type"] == "alleles":
            result["type"] = "custom-allelegraph"
            result["graph"] = convert_allele_vcf(tf.name, None,
                                                 ref_node_padding=params["read_length"],
                                                 ref_node_max_length=params["max_ref_node_length"])
        else:
            result["type"] = "custom-haplotypegraph"
            result["graph"] = convert_haploid_vcf(tf.name,
                                                  None,
                                                  params["reference"],
                                                  ref_node_padding=params["read_length"],
                                                  ref_node_max_length=params["max_ref_node_length"])

        # reformat the result
        if "model_name" in result["graph"]:
            result["ID"] = result["graph"]["model_name"]

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
