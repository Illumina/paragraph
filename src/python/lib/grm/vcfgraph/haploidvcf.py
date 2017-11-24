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

from collections import defaultdict
import re
import pysam
from intervaltree import IntervalTree

from grm.helpers import parse_region


class HaploidVCFGraph(object):
    """
    Class to construct graphs from haploid VCFs. Each sample specifies
    one haplotype which gets labelled as a sequence in the resulting
    paragraph JSON graph.
    """
    # pylint: disable=too-many-instance-attributes
    def __init__(self, chrom="chr"):
        self.ref_split_positions = set()
        self.alts = {}
        self.refs = IntervalTree()
        self.breakends = {}
        self.all_sequences = {"REF"}
        self.chrom = chrom
        self.first_pos = None
        self.last_pos = None

    def dump(self):
        """
        Dump all contents to stdout
        """
        print("Locations: %s:%s-%s" % (str(self.chrom), str(self.first_pos), str(self.last_pos)))
        print("Sequences: %s" % str(list(self.all_sequences)))
        print("Ref-splits: %s" % str(sorted(list(self.ref_split_positions))))
        from pprint import pformat
        print("Alts: %s" % pformat(self.alts))
        print("Reference intervals: %s" % pformat(list(self.refs.items())))
        print("Breakends: %s" % pformat(list(self.breakends)))

    @staticmethod
    def create_from_vcf(ref_file_name,
                        vcf_file_name,
                        chrom=None, start=None, end=None,
                        padding_length=150):
        """
        Create a VCFGraph object from a VCF input file
        :param ref_file_name: Reference Fasta name
        :param vcf_file_name: VCF file name
        :param chrom: chromosome to read from
        :param start: start position for tabix / VCF fetch
        :param end: end position for tabix / VCF fetch
        :param padding_length: length of reference padding
        :return: a VCFGraph
        :rtype: HaploidVCFGraph
        """
        ref_fasta = pysam.FastaFile(ref_file_name)
        vcf = pysam.VariantFile(vcf_file_name)

        graph = HaploidVCFGraph(chrom)

        for sample in vcf.header.samples:  # pylint: disable=no-member
            graph.add_sequence(str(sample))

        first_pos = None
        last_pos = None
        for record in vcf.fetch(chrom, start, end):
            if chrom is None:
                chrom = record.chrom
                graph.chrom = chrom
            elif chrom != record.chrom:
                break

            if first_pos is None:
                first_pos = record.pos

            rec_end = record.stop
            for sample in record.samples:
                for genotype in record.samples[sample]["GT"]:
                    if genotype is not None:
                        if genotype > 0:
                            alt = record.alts[genotype - 1]
                            if alt == "<DEL>":
                                ref_sequence = ref_fasta.fetch(chrom, record.pos - 2, record.stop).upper()
                                graph.add_alt(record.pos - 1, record.stop, ref_sequence, ref_sequence[0], str(sample))
                            elif alt[0] == "<":
                                raise Exception("Unknown symbolic alt: %s:%i - %s" % (chrom, record.pos, alt))
                            elif "SVTYPE" in record.info and record.info["SVTYPE"] == "BND":
                                # We only support forward strand breakends.
                                be_match = re.match(r'([ACGTNXacgtnx]+)([\[\]])([^\[\]]+)([\[\]])', alt)
                                if not be_match:
                                    raise Exception("Unsupported breakend ALT: %s" % alt)
                                pos = record.pos
                                ref_sequence = record.ref
                                ins_sequence = be_match.group(1)
                                be_direction1 = be_match.group(2)
                                be_pos = be_match.group(3)
                                be_direction2 = be_match.group(4)
                                assert be_direction1 == be_direction2
                                if be_direction1 != "[":
                                    raise Exception("Reverse-comp breakends are not supported.")

                                while ref_sequence and ins_sequence and ref_sequence[0] == ins_sequence[0]:
                                    ref_sequence = ref_sequence[1:]
                                    ins_sequence = ins_sequence[1:]
                                    pos += 1

                                be_chrom, be_start, be_end = parse_region(be_pos)
                                assert not be_end
                                if be_chrom != chrom:
                                    raise Exception("Breakends across chromosomes are not supported.")
                                graph.add_breakend(str(sample), pos - 1, ins_sequence, be_start)

                                if not first_pos:
                                    first_pos = record.pos
                                if not last_pos:
                                    last_pos = record.pos

                                first_pos = min([first_pos, be_start, pos])
                                last_pos = max([last_pos, be_start, rec_end])
                            elif re.search(r'[^ACGTNXacgtnx]', alt):
                                raise Exception("Invalid / unsupported ALT: %s" % alt)
                            else:
                                graph.add_alt(record.pos, rec_end, record.ref, alt, str(sample))
                        else:
                            graph.add_ref_support(record.pos, record.stop, str(sample))

            if last_pos is None or last_pos < rec_end:
                last_pos = rec_end

        # add start and end splits
        if first_pos:
            graph.add_ref_split(first_pos - padding_length)
        if last_pos:
            graph.add_ref_split(last_pos + padding_length)

        graph.first_pos = first_pos
        graph.last_pos = last_pos

        # graph.dump()

        return graph

    def add_sequence(self, sequence):
        """ Add a sequence """
        self.all_sequences.add(sequence)

    def add_ref_split(self, pos):
        """ add split for the reference """
        self.ref_split_positions.add(pos)

    def add_alt(self, start, end, ref, alt, sequence):
        """ add alt allele """
        assert len(ref) == end - start + 1

        # trim allele
        while len(alt) > 0 and len(ref) > 0 and ref[-1] == alt[-1]:
            ref = ref[:-1]
            alt = alt[:-1]
            end -= 1
        while len(alt) > 0 and len(ref) > 0 and ref[0] == alt[0]:
            ref = ref[1:]
            alt = alt[1:]
            start += 1

        if not ref and not alt:
            return

        if end < start:
            # insertion
            self.add_ref_split(start)
        else:
            # snp or deletion
            self.add_ref_split(start)
            self.add_ref_split(end + 1)

        key = "%i-%i:%s" % (start, end, alt)
        if key not in self.alts:
            self.alts[key] = {
                "key": key,
                "start": start,
                "end": end,
                "alt": alt,
                "sequences": set(),
                "added-order": len(self.alts)
            }
        if sequence:
            self.alts[key]["sequences"].add(sequence)
            self.add_sequence(sequence)

    def add_ref_support(self, start, end, sequence):
        """
        Tag a piece of reference with a sequence ID
        :param start: start of reference
        :param end: end of reference
        :param sequence: sequence ID
        """
        self.add_ref_split(start)
        self.add_ref_split(end + 1)
        self.refs.addi(start, end + 1, sequence)

    def add_breakend(self, sequence, ref_end, alt="", ref_start=None):
        """ Add simple forward breakend with alt sequence
        :param sequence: sequence name supporting the breakend
        :param ref_end: last base before break
        :param ref_start: reference base to go to
        :param alt: sequence to insert between the breakends
        """
        if not ref_start:
            ref_start = ref_end + 1

        self.add_ref_split(ref_end + 1)
        self.add_ref_split(ref_start)

        be_key = "%i:%i" % (ref_end, ref_start)
        if alt:
            be_key += "-" + alt
        if be_key not in self.breakends:
            self.breakends[be_key] = {
                "end": ref_end,
                "start": ref_start,
                "alt": alt,
                "sequences": set()
            }

        if sequence:
            self.breakends[be_key]["sequences"].add(sequence)
            self.add_sequence(sequence)

    def get_nodes_and_edges(self, ref_paths=False, crossovers=False):
        """ create set of nodes from ref and alt entries
        :param ref_paths: create reference path (otherwise, we'll only create ALT paths)
        :param crossovers: create extra edges for ajacent nodes
        """
        ref_positions = sorted(list(self.ref_split_positions))
        if not ref_positions:
            return [], []
        # need at least 2 reference positions
        assert len(ref_positions) >= 2

        chrom = self.chrom

        nodes_by_name = {}
        ref_starting_at = {}
        ref_ending_at = {}
        sorted_ref_nodes = []
        nodes_starting_at = defaultdict(list)
        nodes_ending_at = defaultdict(list)
        edge_dict = {}

        def edge(nf, nt):
            ekey = nf + "_" + nt
            if ekey not in edge_dict:
                edge_dict[ekey] = {
                    "from": nf,
                    "to": nt,
                    "sequences": set(),
                }
            return edge_dict[ekey]

        # Create alt nodes
        alts_by_sequence = defaultdict(list)
        for alt in sorted(self.alts.values(), key=lambda alt_: (alt_["start"], alt_["added-order"])):
            if alt["alt"]:
                alt["node"] = {"name": alt["key"], "sequence": alt["alt"],
                               "sequences": alt["sequences"],
                               "start": alt["start"],
                               "end": alt["end"]}

                if alt["start"] <= alt["end"]:
                    nodes_ending_at[alt["end"]].append(alt["node"])
                    nodes_starting_at[alt["start"]].append(alt["node"])
                else:
                    nodes_ending_at[alt["start"]].append(alt["node"])
                    nodes_starting_at[alt["end"]].append(alt["node"])

                nodes_by_name[alt["key"]] = alt["node"]

            for s in alt["sequences"]:
                alts_by_sequence[s].append(alt)

        # create ref nodes
        start_pos = ref_positions[0]
        previous_ref_node = None
        for end_pos in ref_positions[1:]:
            assert end_pos > start_pos
            node = {
                "name": "ref-%s:%i-%i" % (chrom, start_pos, end_pos - 1),
                "reference": "%s:%i-%i" % (chrom, start_pos, end_pos - 1),
                "sequences": {"REF"},
                "start": start_pos,
                "end": end_pos - 1
            }

            if ref_paths and previous_ref_node:
                edge(previous_ref_node["name"], node["name"])

            previous_ref_node = node
            nodes_by_name[node["name"]] = node
            ref_starting_at[start_pos] = node
            nodes_starting_at[start_pos].append(node)
            ref_ending_at[end_pos - 1] = node
            nodes_ending_at[end_pos - 1].append(node)
            sorted_ref_nodes.append(node)
            start_pos = end_pos

        if crossovers:
            for pos in range(ref_positions[0], ref_positions[-1] + 1):
                for node1 in nodes_ending_at[pos]:
                    for node2 in nodes_starting_at[pos + 1]:
                        if node1["name"] == node2["name"]:
                            continue
                        edge(node1["name"], node2["name"])

        for s in self.all_sequences:
            alt_tree = IntervalTree()
            alts_this_sequence = alts_by_sequence[s]

            for alt in alts_this_sequence:
                if alt["start"] <= alt["end"]:
                    # end < start => insertion ; these should not interfere with reference nodes
                    alt_tree.addi(alt["start"], alt["end"] + 1, alt)
                    assert alt["start"] - 1 in ref_ending_at
                    assert alt["end"] + 1 in ref_starting_at
                else:
                    # insertion case: check we have the right ref nodes
                    assert alt["end"] in ref_ending_at
                    assert alt["start"] in ref_starting_at

            sorted_ref_this_sequence = [node for node in sorted_ref_nodes
                                        if not alt_tree.overlaps(node["start"], node["end"] + 1)]

            alt_nodes_this_sequence = [alt["node"] for alt in alts_this_sequence
                                       if "node" in alt]
            sorted_nodes_this_sequence = sorted(alt_nodes_this_sequence + sorted_ref_this_sequence,
                                                key=lambda n_: int(n_["start"]))

            previous_node = None
            for current_node in sorted_nodes_this_sequence:
                start_pos = min(current_node["start"], current_node["end"])
                end_pos = max(current_node["start"], current_node["end"])
                ref_supported = len([x for x in self.refs.search(start_pos, end_pos + 1)
                                     if x.data == s]) > 0
                if "sequence" in current_node or ref_supported:
                    current_node["sequences"].add(s)
                if previous_node:
                    e = edge(previous_node["name"], current_node["name"])
                    if (s in previous_node["sequences"] and s in current_node["sequences"]) or previous_node["end"] + 1 < start_pos:
                        e["sequences"].add(s)
                previous_node = current_node

        for be_key, be in self.breakends.items():
            e2 = None
            key = "be_alt_" + be_key
            if be["alt"]:
                alt_sequence = be["alt"]
                nodes_by_name[key] = {"name": key, "sequence": alt_sequence,
                                      "sequences": be["sequences"],
                                      "start": be["start"],
                                      "end": be["end"]}
                e1 = edge(ref_ending_at[be["end"]]["name"], key)
                if be["start"]:
                    e2 = edge(key, ref_starting_at[be["start"]]["name"])
            elif be["start"]:
                e1 = edge(ref_ending_at[be["end"]]["name"], ref_starting_at[be["start"]]["name"])
            else:
                raise Exception("BE %s: must either have ALT or end point." % be_key)

            for s in be["sequences"]:
                e1["sequences"].add(s)
                if e2:
                    e2["sequences"].add(s)

        all_nodes = sorted(nodes_by_name.values(), key=lambda x: min(x["start"], x["end"]))
        for n in all_nodes:
            del n["start"]
            del n["end"]
            if "next_ref" in n:
                del n["next_ref"]
            if n["sequences"] == self.all_sequences:
                del n["sequences"]
            else:
                n["sequences"] = list(n["sequences"])

        all_edges = sorted(edge_dict.values(), key=lambda e_: (e_["from"], e_["to"]))

        for e in all_edges:
            assert e["from"] in nodes_by_name
            assert e["to"] in nodes_by_name
            if e["sequences"] == self.all_sequences or not e["sequences"]:
                del e["sequences"]
            else:
                e["sequences"] = list(e["sequences"])

        return all_nodes, all_edges
