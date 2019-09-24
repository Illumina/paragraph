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
# May 2018
#
# Update VCF file with genotype information from grmpy output
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

from collections import defaultdict
import gzip
import json
import logging
import pysam  # pylint: disable=E0401

from grm.vcfgraph.vcfgraph import VCFGraph


# generate all possible GTs (as described in the VCF SPEC)
def makePLGenotypes(ploidy, alleles, suffix=None, gts=None):
    if not isinstance(gts, list):
        gts = []
    if not isinstance(suffix, list):
        suffix = []
    for allele in range(alleles + 1):
        if ploidy == 1:
            new_suffix = suffix[:]
            new_suffix.insert(0, allele)
            gts.append(new_suffix)
        elif ploidy > 1:
            new_suffix = suffix[:]
            new_suffix.insert(0, allele)
            makePLGenotypes(ploidy - 1, allele, new_suffix, gts)
    return gts


def read_grmpy(grmpyJsonName):
    """
    Read grmpy output
    :param grmpyJsonName: input filename
    :return: grmpy output
    """
    logging.debug("Reading %s", grmpyJsonName)
    if grmpyJsonName.endswith(".gz"):
        f = gzip.open(grmpyJsonName, "rt")
    else:
        f = open(grmpyJsonName, "rt")
    data = json.load(f)

    logging.debug("Ordering by sequences")
    by_sequencename = defaultdict(set)
    by_id = defaultdict(set)
    values = data
    if isinstance(data, dict):
        values = [data]
    for i, d in enumerate(values):
        try:
            sequences = d["graphinfo"]["sequencenames"]
        except KeyError:
            sequences = []

        try:
            ident = d["graphinfo"]["ID"]
        except KeyError:
            ident = None

        if ident:
            by_id[ident].add(i)
        for x in sequences:
            by_sequencename[x].add(i)

    result = {
        "by_id": {k: [values[i] for i in r] for k, r in by_id.items()},
        "by_sequencename": {k: [values[i] for i in r] for k, r in by_sequencename.items()},
    }

    logging.debug("Done reading %s", grmpyJsonName)
    return result


def update_vcf_from_grmpy(inVcfFilename, grmpyOutput, outVcfFilename, sample_names=None):
    inVariantFile = pysam.VariantFile(inVcfFilename)

    # pylint: disable=no-member
    header = inVariantFile.header.copy()

    vcf_samples = list(header.samples)
    if vcf_samples:
        header.add_line('##FORMAT=<ID=OLD_GT,Number=1,Type=String,Description="Previous GT which was replaced by paragraph">')

    if sample_names is None:  # get sample names from vcf
        sample_names = vcf_samples
        if not sample_names:
            raise Exception("Didn't find sample names in either input or VCF. Paragraph cannot output any genotypes!")
    sample_names = set(sample_names)
    sample_names.update(set(vcf_samples))
    added_samples = sample_names.difference(set(vcf_samples))

    sample_names = list(sample_names)
    added_samples = list(added_samples)
    for s in added_samples:
        header.add_sample(s)

    if 'GT' not in header.formats:
        header.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    if 'FT' not in header.formats:
        header.add_line('##FORMAT=<ID=FT,Number=1,Type=String,Description="Filter for genotype">')
    if 'DP' not in header.formats:
        header.add_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total filtered read depth used for genotyping.">')
    if 'AD' not in header.formats:
        header.add_line(
            '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depth for each allele, including the reference.">')
    if 'ADF' not in header.formats:
        header.add_line(
            '##FORMAT=<ID=ADF,Number=R,Type=Integer,Description="Allele depth on forward strand for each allele, including the reference.">')
    if 'ADR' not in header.formats:
        header.add_line(
            '##FORMAT=<ID=ADR,Number=R,Type=Integer,Description="Allele depth on reverse strand for each allele, including the reference.">')
    if 'PL' not in header.formats:
        header.add_line('##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled '
                        'likelihoods for genotypes as defined in the VCF specification">')
    if 'GRMPY_ID' not in header.info:
        header.add_line('##INFO=<ID=GRMPY_ID,Number=1,Type=String,Description="Graph ID '
                        'for linking to genotypes.json.gz; matches record.graphinfo.ID in there.">')
    header.add_line('##FILTER=<ID=BP_DEPTH,Description="One or more breakpoints have abnormal depth">')
    header.add_line('##FILTER=<ID=NO_VALID_GT,Description="No valid genotypes from breakpoints">')
    header.add_line('##FILTER=<ID=CONFLICT,Description="Breakpoints gave different genotypes">')
    header.add_line('##FILTER=<ID=BP_NO_GT,Description="One genotype was missing">')
    header.add_line('##FILTER=<ID=NO_READS,Description="No reads could be retrieved for a breakpoint.">')
    header.add_line(
        '##FILTER=<ID=DEPTH,Description="Poisson depth filter: observed depth deviates too far from Poisson expectation">')
    header.add_line('##FILTER=<ID=UNMATCHED,Description="VCF record could not be matched to a paragraph record.">')
    header.add_line('##FILTER=<ID=MULTIMATCHED,Description="VCF record could not be matched to a paragraph record uniquely.">')
    # pylint: enable=no-member

    bcf_out = pysam.VariantFile(outVcfFilename, 'w', header=header)
    matched = 0
    unmatched = 0
    multimatched = 0
    for raw_record in inVariantFile:
        logging.debug("VCF to graph matching statistics: %d matched, %d unmatched, %d multiple matches.",
                      matched, unmatched, multimatched)
        # reset varIdCounts for every record. This
        # doesn't quite do the right thing when we have
        # multiple VCF variants per block, and where some of these variants
        # start at the same position but come in different records. This could be fixed
        # by collapsing these into single VCF records as alleles before running
        # through paragraph.
        try:
            record = header.new_record(contig=raw_record.chrom, start=raw_record.start, stop=raw_record.stop,
                                       alleles=raw_record.alleles, id=raw_record.id, qual=raw_record.qual, filter=raw_record.filter, info=raw_record.info)
        except:
            raise Exception("Format error in vcf line: " + str(raw_record))

        varIdCounts = defaultdict(int)
        varId = VCFGraph.generate_variant_id(record, varIdCounts)
        alleleIds = [aId for aId, _ in VCFGraph.generate_allele_ids(record, varId)]

        try:
            grmpyRecords = [grmpyOutput["by_id"][raw_record.info["GRMPY_ID"]]]
        except KeyError:
            grmpyRecords = []

        if not grmpyRecords:
            # match by sequence name if we don't have GRMPY_ID
            grmpyRecords = [grmpyOutput["by_sequencename"][aId]
                            for aId in alleleIds if aId in grmpyOutput["by_sequencename"]]

        records = []
        for recordList in grmpyRecords:
            for grmpyRecord in recordList:
                if grmpyRecord not in records:
                    records.append(grmpyRecord)

        if not records:
            record.info["GRMPY_ID"] = "UNMATCHED"
            record.filter.add("UNMATCHED")
            bcf_out.write(record)
            unmatched += 1
            continue
        elif len(records) == 1:
            matched += 1
        else:
            multimatched += 1
            record.info["GRMPY_ID"] = "MULTIPLE:" + ",".join([record["graphinfo"]["ID"] for record in records
                                                              if "graphinfo" in record and "ID" in record["graphinfo"]])
            record.filter.add("MULTIMATCHED")
            bcf_out.write(record)
            continue

        grmpyRecord = records[0]

        try:
            record.info["GRMPY_ID"] = grmpyRecord["graphinfo"]["ID"]
        except KeyError:
            record.info["GRMPY_ID"] = "NOID"

        alleleMap = {"REF": 0, "ALT": 1}
        for ii, aId in enumerate(alleleIds):
            alleleMap[aId] = ii

        num_bpdepth_sample = 0
        for sample in sample_names:
            if vcf_samples:
                if sample in vcf_samples:
                    for k, v in raw_record.samples[sample].items():
                        record.samples[sample][k] = v
                    old_gt = "/".join(sorted([str(val)
                                              if val is not None else "." for val in raw_record.samples[sample]["GT"]]))
                    record.samples[sample]["OLD_GT"] = old_gt
                else:
                    record.samples[sample]["OLD_GT"] = [None]
            record.samples[sample]["GT"] = [None]
            record.samples[sample]["DP"] = None
            record.samples[sample]["FT"] = None
            record.samples[sample]["AD"] = [None] * (1 + len(record.alts))
            record.samples[sample]["ADF"] = [None] * (1 + len(record.alts))
            record.samples[sample]["ADR"] = [None] * (1 + len(record.alts))
            if sample in grmpyRecord["samples"]:
                try:
                    set_record_for_sample(record, sample, grmpyRecord, alleleMap)
                except KeyError:
                    logging.warning("VCF key error for sample " + str(sample) + " at " +
                                    '_'.join(map(str, [record.chrom, record.pos, record.stop])))
                    continue
                if "BP_DEPTH" in record.samples[sample]["FT"] or "BP_NO_GT" in record.samples[sample]["FT"]:
                    num_bpdepth_sample += 1
        if num_bpdepth_sample * 2 > len(grmpyRecord["samples"]):
            record.filter.add("BP_DEPTH")

        bcf_out.write(record)

    logging.info(
        f"VCF to graph matching statistics: {matched} matched, {unmatched} unmatched, {multimatched} multiple matches.")


def set_record_for_sample(record, sample, grmpyRecord, alleleMap):
    """
    set vcf record for the given sample
    """
    gt = grmpyRecord["samples"][sample]["gt"]
    filters = set()
    for f in gt["filters"]:
        filters.add(f)
    gt_to_set = sorted([alleleMap[g] if g in alleleMap else -1 for g in gt["GT"].split("/")])
    gt_to_set = [g if g >= 0 else None for g in gt_to_set]
    if None in gt_to_set:
        # this happens when we cannot match an allele
        filters.add("UNMATCHED")
    else:
        record.samples[sample]["GT"] = gt_to_set
    record.samples[sample]["FT"] = [",".join(list(filters))]
    try:
        record.samples[sample]["DP"] = gt["num_reads"]
    except KeyError:
        record.samples[sample]["DP"] = 0

    ad = grmpyRecord["samples"][sample]["alleles"]
    ads_to_set = [0] * (1 + len(record.alts))
    adfs_to_set = [0] * (1 + len(record.alts))
    adrs_to_set = [0] * (1 + len(record.alts))
    for a in ad.keys():
        ads_to_set[alleleMap[a]] = ad[a]['num_fwd_reads'] + ad[a]['num_rev_reads']
        adfs_to_set[alleleMap[a]] = ad[a]['num_fwd_reads']
        adrs_to_set[alleleMap[a]] = ad[a]['num_rev_reads']
    record.samples[sample]["AD"] = ads_to_set
    record.samples[sample]["ADF"] = adfs_to_set
    record.samples[sample]["ADR"] = adrs_to_set

    ploidy = len(gt_to_set)
    allelecount = len(record.alts)

    gtlist = makePLGenotypes(ploidy, allelecount)
    gtlist_map = {}
    for ii, g in enumerate(gtlist):
        gtlist_map[str(g)] = ii
    pls_to_set = [0] * len(gtlist)
    min_pl = None
    if "GL" not in gt:  # missing genotype
        return
    for name, ll in gt["GL"].items():
        alleles = sorted([alleleMap[a] for a in name.split("/")])
        try:
            phred_l = round(-10 * ll)
        except TypeError:
            phred_l = None
        except OverflowError:
            phred_l = 32768
        phred_l = min(phred_l, 32768)

        if min_pl is None or phred_l < min_pl:
            min_pl = phred_l

        if str(alleles) in gtlist_map:
            pls_to_set[gtlist_map[str(alleles)]] = phred_l

    # normalize, see e.g. https://software.broadinstitute.org/gatk/documentation/article?id=5913
    pls_to_set = [pl - min_pl for pl in pls_to_set]
    record.samples[sample]["PL"] = pls_to_set
