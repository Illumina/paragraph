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
# Class to store reference-based graphs from VCF input (single contig)
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

from collections import namedtuple, defaultdict
import re
import logging
import pysam  # pylint: disable=E0401
import intervaltree

from grm.helpers import parse_region
from grm.vcfgraph.graphContainer import GraphContainer


class NoVCFRecordsException(Exception):
    pass


class VCFGraph:
    """
    Class to construct a more graph-like representation from VCFs.
    Each sample specifies one haplotype (if used)
    Intermediate step used to create paragraph graphs.
    """

    RefInfo = namedtuple("RefInfo", ["haplotypes"])
    AltAllele = namedtuple("AltAllele", ["key", "start", "end", "sequence", "haplotypes"])

    # pylint: disable=too-many-instance-attributes
    def __init__(self, refFile, chrom="chr"):
        self.ref_fasta = pysam.FastaFile(refFile)
        self.alts = {}  # AltAllele.Key -> AltAllele
        self.refs = intervaltree.IntervalTree()
        self.chrom = chrom
        self.first_pos = None
        self.last_pos = None

    def __str__(self):
        """ Summarize graph contents
        """
        from pprint import pformat
        s = [
            "Locations: %s:%s-%s" % (str(self.chrom), str(self.first_pos), str(self.last_pos)),
            "Sequences: %s" % str(list(self.get_haplotypes())),
            "Alt alleles:", pformat(list(self.alts.values())),
            "Reference intervals:", pformat(list(self.get_ref_alleles()))
        ]
        return "\n".join(s)

    @staticmethod
    def generate_variant_id(record, varIdCounts=None):
        """
        Generate a variant ID for a pySAM VCF record
        :param record: a pysam record
        :param varIdCounts: defaultdict dictionary of counts per variant ID:
                            varIdCounts = defaultdict(int)
        :return: variant ID string
        """
        if record.id:
            varId = record.id
            if varIdCounts is not None:
                if varId in varIdCounts:
                    raise Exception(f"Duplicated variant ID: {varId}")
                varIdCounts[varId] = 1
        else:
            varId = f"{record.chrom}:{record.pos}"
            if varIdCounts is not None:
                varIdCounts[varId] += 1
                varId = f"{varId}-{varIdCounts[varId]}"
        return varId

    @staticmethod
    def generate_allele_ids(record, varId):
        return [(f"{varId}:{n}", record.alleles[n]) for n in range(len(record.alleles))]

    @staticmethod
    def create_from_vcf(ref_file_name,
                        vcf_file_name,
                        chrom=None, start=None, end=None,
                        padding_length=150,
                        allele_graph=False):
        """ Create a VCFGraph object from a VCF input file
        :param ref_file_name: Reference Fasta name
        :param vcf_file_name: VCF file name
        :param chrom: chromosome to read from
        :param start: start position for tabix / VCF fetch
        :param end: end position for tabix / VCF fetch
        :param padding_length: length of reference padding
        :param allele_graph: Create an allele graph, rather than haplotype graph
        :return: a VCFGraph
        """
        vcf = pysam.VariantFile(vcf_file_name)
        graph = VCFGraph(ref_file_name, chrom)
        varIdCounts = defaultdict(int)
        record_count = 0
        for record in vcf.fetch(chrom, start, end):
            logging.debug(f"Processing: {str(record).rstrip()}")
            if chrom is None:
                chrom = record.chrom
                graph.chrom = chrom
            elif chrom != record.chrom:
                break
            if graph.first_pos is None:
                graph.first_pos = record.pos
            if graph.last_pos is None or graph.last_pos < record.stop:
                graph.last_pos = record.stop

            varId = VCFGraph.generate_variant_id(record, varIdCounts)
            record_count += 1
            graph.add_record(record, allele_graph, varId)

        if not record_count:
            raise NoVCFRecordsException(f"No VCF records found at {chrom}:{start}-{end}")

        graph.add_ref_support(graph.first_pos - padding_length, graph.last_pos + padding_length)
        # Split read nodes for ALTs to link into (esp. for remote breakends)
        for be in graph.alts.values():
            if graph.first_pos <= be.end <= graph.last_pos:
                graph.refs.slice(be.end + 1)
            else:
                graph.add_ref_support(be.end + 1, be.end + padding_length)
        return graph

    def add_record(self, vcf, allele_graph, varId):
        """ Add one vcfRecord to the graph
        :param vcf: VCF record
        :param allele_graph: Use all alleles (from individual VCF entries), rather than haplotypes (from VCF samples)
        :param varId: Globally unique ID to use for this variant (e.g. vcf.id if not None)
        """
        if not allele_graph:
            try:
                samples = {s.name: s.alleles[0] for s in vcf.samples.values() if None not in s.alleles}
            except:  # pylint: disable=bare-except
                # when no samples are present, don't label the haplotype, use allele labels instead
                samples = {x: y for x, y in VCFGraph.generate_allele_ids(vcf, varId)}
        else:
            samples = {x: y for x, y in VCFGraph.generate_allele_ids(vcf, varId)}

        refSamples = set(s for s in samples if samples[s] == vcf.ref)
        self.add_ref_support(vcf.pos, vcf.stop, refSamples, vcf.alleles)
        vcfalts = vcf.alts if vcf.alts else []
        for alt in vcfalts:
            alt_samples = set(s for s in samples if samples[s] == alt)
            if not (alt_samples or allele_graph):
                continue  # Skip alleles not used in any haplotype
            if alt == "<DEL>":
                ref_sequence = self.ref_fasta.fetch(self.chrom, vcf.pos-1, vcf.stop).upper()
                self.add_alt(vcf.pos, vcf.stop, ref_sequence, "", alt_samples)
            elif alt == "<INS>":
                i_ref_start = vcf.pos
                i_ref_end = vcf.stop
                ref_sequence = self.ref_fasta.fetch(self.chrom, i_ref_start-1, i_ref_end).upper()
                if len(ref_sequence) != i_ref_end - i_ref_start + 1:
                    raise Exception(f"Cannot retrieve reference sequence for {self.chrom}:{i_ref_start}-{i_ref_end}; "
                                    "are you using the correct FASTA reference?")
                if "SEQ" not in vcf.info:
                    raise Exception(f"No insertion sequence specified for <INS> at {self.chrom}:{i_ref_start}-{i_ref_end}; ")
                self.add_alt(i_ref_start, i_ref_end, ref_sequence, ref_sequence[0] + vcf.info["SEQ"],
                             alt_samples, refSamples)
            elif alt[0] == "<":
                raise Exception("Unknown symbolic alt: %s:%i - %s" % (self.chrom, vcf.pos, alt))
            elif "SVTYPE" in vcf.info and vcf.info["SVTYPE"] == "BND":
                ins_sequence, be_start = self._parse_breakend(alt)
                self.add_breakend(vcf.pos, vcf.ref, be_start, alt_samples, ins_sequence, refSamples)
            elif re.search(r'[^ACGTNXacgtnx]', alt):
                raise Exception("Invalid / unsupported ALT: %s" % alt)
            else:
                self.add_alt(vcf.pos, vcf.stop, vcf.ref, alt, alt_samples, refSamples)

    def add_ref_support(self, start, end, haplos=(), alleles=None):
        """ Tag a piece of reference with a haplotype
        Don't apply haplotype label to leading padding sequence shared by all alleles
        :param start: start of reference
        :param end: end of reference (inclusive)
        :param haplos: Haplotypes supporting reference
        :param alleles: Other alleles in this VCF entry
        """
        pad = 0
        if alleles:
            minLen = min(len(a) for a in alleles)
            while pad < minLen and all(alleles[0][pad] == a[pad] for a in alleles):
                pad += 1
        assert start + pad <= end + 1
        logging.debug(f"Adding REF: {start}-{end} Haplotypes:{haplos}")

        if pad > 0:
            # Create full-length reference block, but only apply haplo label to non-padding bases
            self.refs.addi(start, end + 1, VCFGraph.RefInfo(set()))
            logging.debug(f"Skipping {pad} ref-padding bases")
            if haplos and start + pad <= end:
                self.refs.addi(start + pad, end + 1, VCFGraph.RefInfo(set(haplos)))
        else:
            self.refs.addi(start, end + 1, VCFGraph.RefInfo(set(haplos)))

    def get_ref_alleles(self):
        """ Split reference intervals into non-overlapping pieces, preserving haplotypes
        :return List of reference intervals (sorted by start/stop)
        """
        self.refs.split_overlaps()
        lastNode = None
        for ref in sorted(self.refs):
            if not lastNode:
                lastNode = ref
            elif not ref.range_matches(lastNode):
                yield lastNode
                lastNode = ref
            else:
                haplos = lastNode.data.haplotypes.union(ref.data.haplotypes)
                lastNode = lastNode._replace(data=VCFGraph.RefInfo(haplos))
        if lastNode:
            yield lastNode

    def add_alt(self, start, end, ref, alt, haplos=(), other_haplos=()):
        """ Add one alt allele to the graph
        Trim of leading padding bases (shared with ref)
        For insertions add a 'bypass allele' parallel to the insertion for haplotypes
        typed for a different allele at this locus (e.g. ref)
        :param start: Start coordinate of allele in reference
        :param end: End coordiante of allele in reference (inclusive)
        :param ref: Reference allele sequence
        :param alt: Alt allele sequence
        :param haplos: Haplotypes with alt allele
        :param other_haplos: Haplotypes typed for another allele at this locus
        """
        assert len(ref) == end - start + 1
        # trim alt allele
        alt_start, alt_end = start, end
        while alt and ref and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            alt_start += 1
        if alt_start > start:
            # Don't count padding bases as a 'reference call', i.e. don't label with haplotype
            self.add_ref_support(start, alt_start - 1)
        while alt and ref and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]
            alt_end -= 1
        assert alt_end > 0
        if alt_start <= alt_end < end:
            # add reference support when we have trimmed unless it's an insertion
            self.add_ref_support(alt_end + 1, end, haplos)
        assert ref or alt

        self._addAlt(alt_start, alt_end, alt, haplos)
        # Add 'bypass' node for insertions to prevent other_haplos from using insertion
        # This becomes an edge bypassing the insertion labeled with other_haplos later
        if other_haplos and alt_start > alt_end:
            self._addAlt(alt_start, alt_end, "", other_haplos)

    def _parse_breakend(self, alt):
        """ Parse remote breakend info from BND VCF record
        :param record: VCF record
        :param alt: Alt allele from VCF record (BND)
        :return Inserted sequence and position of first base after remote breakend
        """
        # We only support forward strand breakends.
        be_match = re.match(r'([ACGTNXacgtnx]+)([\[\]])([^\[\]]+)([\[\]])', alt)
        if not be_match:
            raise Exception("Unsupported breakend ALT: %s" % alt)
        ins_sequence = be_match.group(1)
        be_direction1 = be_match.group(2)
        be_pos = be_match.group(3)
        be_direction2 = be_match.group(4)
        be_chrom, be_start, be_end = parse_region(be_pos)
        if be_direction1 != "[" or be_direction2 != "[":
            raise Exception("Reverse-comp breakends are not supported.")
        assert not be_end
        if be_chrom != self.chrom:
            raise Exception("Breakends across chromosomes are not supported.")
        return ins_sequence, be_start

    def add_breakend(self, pos, ref_seq, end, haplos=(), ins_seq="", ref_haplos=()):
        """ Add simple breakend allele to the graph
        :param pos: Base of breakend
        :param ref_seq: Sequence of ref allele
        :param end: Coordinate of first base after remote breakend (same chrom)
        :param haplos: Haplotypes with alt allele
        :param ins_seq: Inserted breakpoint sequence
        :param ref_haplos: Haplotypes with reference allele
        """
        # trim alt allele
        alt_start = pos
        while ref_seq and ins_seq and ref_seq[0] == ins_seq[0]:
            ref_seq = ref_seq[1:]
            ins_seq = ins_seq[1:]
            alt_start += 1
        assert alt_start != end - 1
        # Reference node covering the entire skipped sequence not part of any allele. Will likely be split later.
        self.add_ref_support(pos, end - 1)
        self._addAlt(alt_start, end - 1, ins_seq, haplos)
        # 'Pseudo-insertions' that create edges spanning the breakends labeled with reference sequences.
        self._addAlt(alt_start, alt_start - 1, "", ref_haplos)
        self._addAlt(end, end - 1, "", ref_haplos)

    def _addAlt(self, start, end, seq, haplos=()):
        key = f"{start}-{end}:{seq}"
        logging.debug(f"Adding ALT: {key} Seqs:{haplos}")
        if key not in self.alts:
            self.alts[key] = VCFGraph.AltAllele(key, start, end, seq, set())
        self.alts[key].haplotypes.update(haplos)


    def get_haplotypes(self):
        """
        :return All haplotype IDs (i.e. vcf sample names) present in graph
        """
        return (
            set(s for a in self.alts.values() for s in a.haplotypes) |
            set(s for i in self.refs for s in i.data.haplotypes)
        ).difference([None])

    def get_graph(self, allele_graph=False):
        """ Create the paragraph representation of nodes and edges for this graph
        :param alleleGraph: create edges between any compatible allele pair (rather
                            than just following reference and given haplotypes)
        :return GraphContainer object
        """
        logging.info("Creating output graph")
        graph = GraphContainer()
        # create ref nodes
        pnode = None
        for ref in self.get_ref_alleles():
            node = graph.add_refNode(self.chrom, ref.begin, ref.end - 1, ref.data.haplotypes)
            if pnode:
                assert pnode["end"] + 1 == node["start"]
                graph.add_edge(pnode, node)
            pnode = node
        # Create alt nodes
        for alt in self.alts.values():
            graph.add_altNode(self.chrom, alt.start, alt.end, alt.sequence, alt.haplotypes)

        # Create edges connecting nodes along a haplotype (or allele in alleleGraph mode)
        for haplo in self.get_haplotypes():
            nodes = graph.nodes_by_haplo(haplo)
            logging.info(f"Linking nodes in sequence {haplo}:\t{', '.join(n['name'] for n in nodes)}")
            pnode = None
            for node in nodes:
                if pnode:
                    if pnode["end"] == node["start"] - 1:
                        graph.add_edge(pnode, node, [haplo])
                    pnode_is_ref_dummy = pnode["end"] == pnode["start"] - 1 and not pnode["sequence"]
                    pnode_ends_before_node = pnode["end"] < node["start"] and pnode["start"] < node["start"]
                    if not pnode_is_ref_dummy and not pnode_ends_before_node:
                        raise Exception(f"Inconsistent nodes for haplotype {haplo}: {pnode['name']}, {node['name']}")
                pnode = node

        # In alleleGraph mode link each alt node to all neighboring nodes
        # In haplotype mode link nodes without in/out edges to reference
        for node in graph.altNodes():
            if allele_graph or not any(graph.inEdges(node)):
                graph.add_edge(graph.refNode_ending_at[node["chrom"], node["start"] - 1], node)
            if not any(graph.outEdges(node)):
                graph.add_edge(node, graph.refNode_starting_at[node["chrom"], node["end"] + 1])
            if allele_graph:
                isInsertion = node["end"] < node["start"]
                for n in graph.nodes_starting_at[node["end"]+1]:
                    # Don't loop by connecting multiple insertions at the same position
                    if not (isInsertion and n["end"] < n["start"]):
                        graph.add_edge(node, n)

        # For nodes that do not have determined in/out edges for a given haplotype
        # label all in/out edges as compatible with that haplotype
        # excluding edges that connect to another allele at the same vcfVariant (e.g. insertions)
        for haplo in self.get_haplotypes():
            for node in graph.nodes_by_haplo(haplo):
                if not any(graph.inEdges(node, haplo)):
                    for e in graph.inEdges(node):
                        graph.add_edge(graph.nodes[e["from"]], node, [haplo])
                assert any(graph.inEdges(node, haplo))
                if not any(graph.outEdges(node, haplo)):
                    for e in graph.outEdges(node):
                        graph.add_edge(node, graph.nodes[e["to"]], [haplo])
        return graph
