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
import os
import argparse
import subprocess
import tempfile
import shlex
import pysam

import findgrm  # pylint: disable=unused-import

from grm.msa import *  # pylint: disable=wildcard-import,unused-wildcard-import

VCFHEADER = """##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=END,Number=1,Type=Integer,Description="End of SV/homref">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">"""


def main():
    parser = argparse.ArgumentParser("vcf2paragraph.py")

    parser.add_argument("input", help="Input Fasta MSA file", nargs=1)
    parser.add_argument("output", help="Output VCF file", nargs=1)

    parser.add_argument("--reference-chr", dest="reference_chr",
                        help="Reference chromosome", default="chr")
    parser.add_argument("--reference-start", help="Reference start pos",
                        type=int, default=0, dest="reference_start")
    parser.add_argument("--reference-sequence", help="Reference FASTA for checking REF",
                        type=str, dest="ref", required=True)
    parser.add_argument(
        "--bcftools", help="Path to bcftools", default="bcftools")

    args = parser.parse_args()

    fasta_seqs = read_msa_fasta(args.input[0])

    s0 = list(fasta_seqs.values())[0]

    ref_fasta = pysam.FastaFile(args.ref)
    s0_nodels = s0.replace("-", "").upper()
    full_reference_seq = ref_fasta.fetch(args.reference_chr, args.reference_start - 1,
                                         args.reference_start + len(s0_nodels) - 1).upper()
    assert s0_nodels == full_reference_seq

    args.reference_end = args.reference_start + len(s0_nodels) - 1

    tf = tempfile.NamedTemporaryFile(suffix=".vcf", delete=False, mode="wt")
    print(VCFHEADER, file=tf)
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s" % list(fasta_seqs.items())[0][0], file=tf)
    if args.reference_start and args.reference_end > args.reference_start:
        print("%s\t%i\t.\t%s\t.\t.\t.\tEND=%i\tGT\t0" % (args.reference_chr,
                                                         args.reference_start,
                                                         full_reference_seq[0],
                                                         args.reference_end),
              file=tf)
    tf.close()
    subprocess.check_call(
        " ".join([args.bcftools, "view", shlex.quote(tf.name), "-o", shlex.quote(tf.name) + ".gz", "-O", "z"]),
        shell=True)
    os.remove(tf.name)
    subprocess.check_call(" ".join([args.bcftools, "index", shlex.quote(tf.name) + ".gz"]), shell=True)

    to_merge = [tf.name + ".gz"]

    for n, s in list(fasta_seqs.items())[1:]:
        pv = pairwise_variants(s0, s, args.reference_start - 1)
        tf = tempfile.NamedTemporaryFile(suffix=".vcf", delete=False, mode="wt")
        print(VCFHEADER, file=tf)
        print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s" % n, file=tf)
        i = 0
        for variant in sorted(pv, key=lambda x: x["start"]):
            i += 1
            if ref_fasta:
                ref_allele = ref_fasta.fetch(args.reference_chr, variant["start"],
                                             variant["start"] + len(variant["ref"])).upper()
                if ref_allele != variant["ref"]:
                    ref_allele = ref_fasta.fetch(args.reference_chr, variant["start"] - 10,
                                                 variant["start"] - 2).lower() + \
                                 ref_allele + ref_fasta.fetch(args.reference_chr,
                                                              variant["start"] + len(variant["ref"]),
                                                              variant["start"] + len(variant["ref"]) + 10).lower()
                    print("%i: %s:%i: ref_allele %s != variant[\"ref\"] %s" % (i, args.reference_chr,
                                                                               variant["start"] + 1,
                                                                               ref_allele,
                                                                               variant["ref"]))
            print("%s\t%i\t.\t%s\t%s\t.\tPASS\t.\tGT\t1" % (
                args.reference_chr, variant["start"] + 1, variant["ref"], variant["alt"]), file=tf)
        tf.close()
        subprocess.check_call(
            " ".join([args.bcftools, "view", shlex.quote(tf.name), "-o", shlex.quote(tf.name) + ".gz", "-O", "z"]),
            shell=True)
        os.remove(tf.name)
        subprocess.check_call(" ".join([args.bcftools, "index", shlex.quote(tf.name) + ".gz"]), shell=True)
        to_merge.append(tf.name + ".gz")

    subprocess.check_call(" ".join([args.bcftools, "merge" if len(to_merge) > 1 else "view",
                                    "-o", args.output[0]] +
                                   list(map(shlex.quote, to_merge))), shell=True)
    for f in to_merge:
        os.remove(f)


if __name__ == "__main__":
    main()
