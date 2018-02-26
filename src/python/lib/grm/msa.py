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
# Functions for dealing with multiple sequence alignments
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

from collections import OrderedDict
from pprint import pprint


def read_msa_fasta(filename):
    """ Read a MSA FASTA file
    :return: dictionary with sequences
    :rtype: OrderedDict
    """

    sequences = OrderedDict()

    with open(filename) as ff:
        current = {}
        for i, line in enumerate(ff):
            line = line.strip()
            if line[0] == ">":
                if current and "name" in current:
                    sequences[current["name"]] = current["sequence"]
                    current = {}
                current["name"], _, current["desc"] = line[1:].partition(" ")
                current["sequence"] = ""
            elif current and "name" in current:
                current["sequence"] += line
            elif line:
                raise Exception("Error in line %i : no fasta sequence name" % i)

    if current and "name" in current:
        sequences[current["name"]] = current["sequence"]

    last = None
    for s in sequences.values():
        if last is None:
            last = len(s)
        assert len(s) == last

    return sequences


def pairwise_variants(ref, alt, offset=0):
    """ Get pairwise variant calls for aligned REF and ALT sequences """

    assert len(ref) == len(alt)

    variants = []
    ref_start = -1
    current_ref_pos = 0
    deleted_sequence = ""
    inserted_sequence = ""
    for i in range(len(ref)):  # pylint: disable=C0200
        if alt[i] != "-" and ref[i] != '-' and ref_start >= 0 and (inserted_sequence or deleted_sequence):
            variants.append({"start": ref_start + offset, "end": ref_start + offset + len(deleted_sequence) - 1,
                             "ref": deleted_sequence, "alt": inserted_sequence})
            inserted_sequence = ""
            deleted_sequence = ""
            ref_start = current_ref_pos

        if ref[i] == alt[i]:
            if ref[i] != "-":
                current_ref_pos += 1
            continue

        if ref[i] == "-":
            inserted_sequence += alt[i]
        elif alt[i] == "-":
            deleted_sequence += ref[i]
            current_ref_pos += 1
        else:
            deleted_sequence += ref[i]
            inserted_sequence += alt[i]
            current_ref_pos += 1

        if len(deleted_sequence) == 1:
            ref_start = current_ref_pos - 1

    if ref_start >= 0 and (inserted_sequence or deleted_sequence):
        variants.append({"start": ref_start + offset, "end": ref_start + len(ref) - 1 + offset,
                         "ref": deleted_sequence, "alt": inserted_sequence})

    variants = sorted(variants, key=lambda v: v["start"])
    last_end = offset

    ref = ref.replace("-", "")

    for variant in variants:
        assert variant["ref"] != variant["alt"]

        def pad(variant_to_pad):
            if variant_to_pad["start"] >= last_end and not variant_to_pad["alt"]:
                variant_to_pad["start"] -= 1
                variant_to_pad["ref"] = ref[variant_to_pad["start"] - offset] + variant_to_pad["ref"]
                variant_to_pad["alt"] = variant_to_pad["ref"][0] + variant_to_pad["alt"]
            elif variant_to_pad["start"] >= last_end and not variant_to_pad["ref"]:
                variant_to_pad["start"] -= 1
                variant_to_pad["ref"] = ref[variant_to_pad["start"] - offset]
                variant_to_pad["alt"] = variant_to_pad["ref"] + variant_to_pad["alt"]

        pad(variant)
        while variant["start"] > last_end and variant["ref"] and variant["alt"] and variant["ref"][-1] == variant["alt"][-1]:
            variant["end"] -= 1
            variant["ref"] = variant["ref"][:-1]
            variant["alt"] = variant["alt"][:-1]
            pad(variant)
        pad(variant)
        last_end = variant["end"]

    pprint(variants)
    return variants
