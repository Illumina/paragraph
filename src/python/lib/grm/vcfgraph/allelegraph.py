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

from grm.helpers import parse_region


class AlleleGraph(object):
    """
    Class to construct graphs from VCFs using allele information.
    Alleles are added as alternate branches to the graph
    """
    # pylint: disable=too-many-instance-attributes
    def __init__(self, chrom="chr"):
        self.ref_split_positions = set()
        self.alts = {}
        self.breakends = {}
        self.all_sequences = {}
        self.chrom = chrom
        self.first_pos = None
        self.last_pos = None
        self.vcflist = []

    def dump(self):
        """
        Dump all contents to stdout
        """
        print("Locations: %s:%s-%s" % (str(self.chrom), str(self.first_pos), str(self.last_pos)))
        print("Sequences: %s" % str(list(self.all_sequences)))
        print("Ref-splits: %s" % str(sorted(list(self.ref_split_positions))))
        from pprint import pformat
        print("Alts: %s" % pformat(self.alts))
        print("Breakends: %s" % pformat(list(self.breakends)))

    @staticmethod
    def create_from_vcf(vcf_file_name,
                        chrom=None, start=None, end=None,
                        padding_length=150):
        """
        Create a graph from a VCF input file
        :param vcf_file_name: VCF file name
        :param chrom: chromosome to read from
        :param start: start position for tabix / VCF fetch
        :param end: end position for tabix / VCF fetch
        :param padding_length: length of reference padding
        :return: a VCFGraph
        :rtype: AlleleGraph
        """
        vcf = pysam.VariantFile(vcf_file_name)

        graph = AlleleGraph(chrom)

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
            vcfinfo = {
                "chrom": chrom,
                "pos": record.pos,
                "end": rec_end,
                "ref": record.ref,
                "alts": [],
                "id": record.id,
            }
            for alt in record.alts:
                vcfinfo["alts"].append(alt)
                if alt == "<DEL>":
                    graph.add_breakend(record.pos - 1, "", rec_end + 1)
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
                    graph.add_breakend(pos - 1, ins_sequence, be_start)

                    if not first_pos:
                        first_pos = record.pos
                    if not last_pos:
                        last_pos = record.pos

                    first_pos = min([first_pos, be_start, pos])
                    last_pos = max([last_pos, be_start, rec_end])
                elif re.search(r'[^ACGTNXacgtnx]', alt):
                    raise Exception("Invalid / unsupported ALT: %s" % alt)
                else:
                    graph.add_alt(record.pos, rec_end, record.ref, alt)
            graph.vcflist.append(vcfinfo)
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

    def add_ref_split(self, pos):
        """ add split for the reference """
        self.ref_split_positions.add(pos)

    def add_alt(self, start, end, ref, alt):
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
            self.add_ref_split(end)
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
                "sequences": [],
                "added-order": len(self.alts)
            }

    def add_breakend(self, ref_end, alt="", ref_start=None):
        """ Add simple forward breakend with alt sequence
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
                "sequences": [],
            }

    def get_nodes_and_edges(self):
        """ create set of nodes from ref and alt entries
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
                    "sequences": [],
                }
            return edge_dict[ekey]

        # create ref nodes
        start_pos = ref_positions[0]
        previous_ref_node = None

        for i, end_pos in enumerate(ref_positions[1:]):
            assert end_pos > start_pos
            node = {
                "name": "REF-%i-%i-%i" % (i + 1, start_pos, end_pos),
                "reference": "%s:%i-%i" % (chrom, start_pos, end_pos - 1),
                "sequences": ["REF"],
                "start": start_pos,
                "end": end_pos - 1
            }

            if previous_ref_node:
                edge(previous_ref_node["name"], node["name"])["sequences"] = ["REF"]

            previous_ref_node = node
            nodes_by_name[node["name"]] = node
            ref_starting_at[start_pos] = node
            nodes_starting_at[start_pos].append(node)
            ref_ending_at[end_pos - 1] = node
            nodes_ending_at[end_pos - 1].append(node)
            sorted_ref_nodes.append(node)
            start_pos = end_pos

        # Create alt nodes
        for i, alt in enumerate(sorted(self.alts.values(), key=lambda alt_: (alt_["start"], alt_["added-order"]))):
            if alt["alt"]:
                alt["node"] = {"name": "ALT-%i-%i-%i" % (i + 1, alt["start"], alt["end"]),
                               "sequence": alt["alt"],
                               "sequences": alt["sequences"],
                               "start": alt["start"],
                               "end": alt["end"]}

                if alt["start"] <= alt["end"]:
                    nodes_ending_at[alt["end"]].append(alt["node"])
                    nodes_starting_at[alt["start"]].append(alt["node"])
                else:
                    nodes_ending_at[alt["end"]].append(alt["node"])
                    nodes_starting_at[alt["end"] + 1].append(alt["node"])

                nodes_by_name[alt["node"]["name"]] = alt["node"]
            else:
                for node1 in nodes_ending_at[alt["start"] - 1]:
                    for node2 in nodes_starting_at[alt["end"] + 1]:
                        edge(node1["name"], node2["name"])

        for pos in range(ref_positions[0], ref_positions[-1] + 1):
            for node1 in nodes_ending_at[pos]:
                for node2 in nodes_starting_at[pos + 1]:
                    if node1["name"] == node2["name"]:
                        continue
                    edge(node1["name"], node2["name"])

        be_num = 0
        for be_key, be in self.breakends.items():
            e2 = None
            key = "BE-%i-%i-%i" % (be_num, be["start"], be["end"])
            if be["alt"]:
                be_num += 1
                alt_sequence = be["alt"]
                nodes_by_name[key] = {"name": key,
                                      "sequence": alt_sequence,
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

        all_edges = sorted(edge_dict.values(), key=lambda e_: (e_["from"], e_["to"]))

        fadj = defaultdict(list)
        for e in all_edges:
            assert e["from"] in nodes_by_name
            assert e["to"] in nodes_by_name
            fadj[e["from"]].append(e["to"])

        self.all_sequences = ["REF"]

        def edges_without_sequences():
            count = 0
            for iedge in edge_dict.values():
                if not iedge["sequences"]:
                    count += 1
            return count

        s_id = 1
        while edges_without_sequences():
            this_node = all_nodes[0]
            s_name = "S%i" % s_id
            this_node["sequences"].append(s_name)
            self.all_sequences.append(s_name)
            while fadj[this_node["name"]]:
                this_node_name = this_node["name"]
                next_node_name = fadj[this_node_name][0]
                edge_n_sequences = len(edge(this_node_name, next_node_name)["sequences"])
                for next_candidate in fadj[this_node["name"]]:
                    if len(edge(this_node_name, next_candidate)["sequences"]) < edge_n_sequences:
                        next_node_name = next_candidate
                        edge_n_sequences = len(edge(this_node_name, next_candidate)["sequences"])
                next_node = nodes_by_name[next_node_name]
                next_node["sequences"].append(s_name)
                edge(this_node["name"], next_node_name)["sequences"].append(s_name)
                this_node = next_node
            s_id += 1

        all_edges = sorted(edge_dict.values(), key=lambda e_: (e_["from"], e_["to"]))

        return all_nodes, all_edges
