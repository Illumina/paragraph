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

import os
import sys
import json
import unittest
import tempfile
import gzip
import difflib

GRMPY_ROOT = os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", ".."))

sys.path.append(os.path.join(GRMPY_ROOT,
                             "src", "python", "bin"))

GRMPY_INSTALL = None
if "GRMPY_INSTALL" in os.environ:
    GRMPY_INSTALL = os.environ["GRMPY_INSTALL"]

hg38_locations = [
    "/Users/pkrusche/workspace/human_genome/hg38.fa",
    "/illumina/sync/igenomes/Homo_sapiens/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa",
]
if "HG38" in os.environ and os.environ["HG38"]:
    hg38_locations.insert(0, os.environ["HG38"])
HG38 = None
for f in hg38_locations:
    if os.path.exists(f) and os.path.isfile(f):
        HG38 = f
        break


class ADict(dict):
    def __init__(self, *args, **kwargs):
        super(ADict, self).__init__(*args, **kwargs)
        self.__dict__ = self


class TestMultiGrmpy(unittest.TestCase):
    def setUp(self):
        self.input_json = os.path.join(GRMPY_ROOT, "share", "test-data", "round-trip-genotyping", "candidates.json")
        self.input_vcf = os.path.join(GRMPY_ROOT, "share", "test-data", "round-trip-genotyping", "candidates.vcf")
        self.reference = os.path.join(GRMPY_ROOT, "share", "test-data", "round-trip-genotyping", "dummy.fa")
        self.manifest = os.path.join(GRMPY_ROOT, "share", "test-data", "round-trip-genotyping", "samples.txt")

        self.swaps_input_vcf = os.path.join(GRMPY_ROOT, "share", "test-data", "genotyping_test_2", "swaps.vcf")
        self.swaps_expected_genotypes_json = os.path.join(GRMPY_ROOT, "share", "test-data", "genotyping_test_2",
                                                          "expected-genotypes.json")
        self.swaps_expected_genotypes_vcf = os.path.join(GRMPY_ROOT, "share", "test-data", "genotyping_test_2",
                                                         "expected-genotypes.vcf")
        self.swaps_expected_variants_json = os.path.join(GRMPY_ROOT, "share", "test-data", "genotyping_test_2",
                                                         "expected-variants.json")
        self.swaps_reference = os.path.join(GRMPY_ROOT, "share", "test-data", "genotyping_test_2", "swaps.fa")
        self.swaps_manifest = os.path.join(GRMPY_ROOT, "share", "test-data", "genotyping_test_2", "samples.txt")

        self.hg38_input_vcf = os.path.join(GRMPY_ROOT, "share", "test-data", "paragraph", "pg-het-ins",
                                           "pg-het-ins.vcf")
        self.hg38_expected_genotypes_json = os.path.join(GRMPY_ROOT, "share", "test-data", "paragraph", "pg-het-ins",
                                                         "genotypes_expected.json")
        self.hg38_expected_variants_json = os.path.join(GRMPY_ROOT, "share", "test-data", "paragraph", "pg-het-ins",
                                                        "variants_expected.json")
        self.hg38_reference = HG38
        self.hg38_manifest = os.path.join(GRMPY_ROOT, "share", "test-data", "paragraph", "pg-het-ins", "manifest.txt")

    @unittest.skipIf(not GRMPY_INSTALL,
                     "No compiled/install path specified, please set the GRMPY_INSTALL variable to point to an installation.")
    def test_multigrmpy(self):
        import multigrmpy

        with tempfile.TemporaryDirectory() as output_dir:
            args = ADict({
                "input": self.input_vcf,
                "manifest": self.manifest,
                "reference": self.reference,
                "output": output_dir,
                "grmpy": os.path.join(os.path.dirname(__file__), "module-wrapper.sh") + " " + os.path.join(
                    GRMPY_INSTALL, "bin", "grmpy"),
                "verbose": False,
                "quiet": True,
                "logfile": None,
                "infer_read_haplotypes": False,
                "write_alignments": False,
                "graph_sequence_matching": True,
                "klib_sequence_matching": False,
                "kmer_sequence_matching": False,
                "bad_align_uniq_kmer_len": 0,
                "threads": 1,
                "sample_threads": 1,
                "genotyping_parameters": None,
                "max_reads_per_event": 10000,
                "split_type": "lines",
                "read_length": 150,
                "max_ref_node_length": 1000,
                "ins_info_key": "SEQ",
                "graph_type": "alleles",
                "alt_splitting": True,
                "retrieve_reference_sequence": False,
                "scratch_dir": None,
                "keep_scratch": True,
            })
            output_json_path = os.path.join(output_dir, "genotypes.json.gz")
            multigrmpy.run(args)
            with gzip.open(output_json_path, 'rt') as result_json:
                observed = json.load(result_json)
                for item in observed:
                    if item["graphinfo"]["ID"] == "test-ins":
                        self.assertEqual(item["samples"]["sample1"]["gt"]["GT"], "S1/S1")
                        self.assertEqual(item["samples"]["sample2"]["gt"]["GT"], "./.")
                    if item["graphinfo"]["ID"] == "test-del":
                        self.assertEqual(item["samples"]["sample1"]["gt"]["GT"], "./.")
                        self.assertEqual(item["samples"]["sample2"]["gt"]["GT"], "S1/S1")

    @unittest.skipIf(not GRMPY_INSTALL,
                     "No compiled/install path specified, please set the GRMPY_INSTALL variable to point to an installation.")
    def test_multigrmpy_expected_genotypes(self):
        import multigrmpy

        with tempfile.TemporaryDirectory() as output_dir:
            args = ADict({
                "input": self.swaps_input_vcf,
                "manifest": self.swaps_manifest,
                "reference": self.swaps_reference,
                "output": output_dir,
                "grmpy": os.path.join(os.path.dirname(__file__), "module-wrapper.sh") + " " + os.path.join(
                    GRMPY_INSTALL, "bin", "grmpy"),
                "verbose": False,
                "quiet": True,
                "logfile": None,
                "infer_read_haplotypes": False,
                "write_alignments": False,
                "graph_sequence_matching": True,
                "klib_sequence_matching": False,
                "kmer_sequence_matching": False,
                "bad_align_uniq_kmer_len": 0,
                "threads": 1,
                "sample_threads": 1,
                "genotyping_parameters": None,
                "max_reads_per_event": 10000,
                "split_type": "lines",
                "read_length": 150,
                "max_ref_node_length": 1000,
                "ins_info_key": "SEQ",
                "graph_type": "alleles",
                "alt_splitting": True,
                "retrieve_reference_sequence": False,
                "scratch_dir": None,
                "keep_scratch": False,
            })
            multigrmpy.run(args)

            output_json_path = os.path.join(output_dir, "genotypes.json.gz")
            with gzip.open(output_json_path, 'rt') as result_json:
                observed = json.load(result_json)

            with open(self.swaps_expected_genotypes_json, 'rt') as expected_json:
                expected = json.load(expected_json)

            # event ordering is not guaranteed
            observed = sorted(observed, key=lambda x: x["graphinfo"]["ID"])
            expected = sorted(expected, key=lambda x: x["graphinfo"]["ID"])

            match = True
            for i, o in enumerate(observed):
                if o["samples"]["SWAPS"]["gt"]["GT"] != expected[i]["samples"]["SWAPS"]["gt"]["GT"]:
                    match = False
                    break

            if not match:
                observed_lines = json.dumps(observed, sort_keys=True, indent=4, separators=(',', ': '))
                expected_lines = json.dumps(expected, sort_keys=True, indent=4, separators=(',', ': '))
                for line in difflib.context_diff(expected_lines.split("\n"), observed_lines.split("\n"),
                                                 fromfile=self.swaps_expected_genotypes_json, tofile=output_json_path):
                    sys.stderr.write(line + "\n")

                with open("test_swaps_genotypes.json", "wt") as out_file:
                    json.dump(observed, out_file, sort_keys=True, indent=4, separators=(',', ': '))

                raise Exception("Swaps test genotyping output doesn't match! If this is expected and new behavior, "
                                "cp test_swaps_genotypes.json %s" % self.swaps_expected_genotypes_json)

            output_vcf_path = os.path.join(output_dir, "genotypes.vcf.gz")

            with gzip.open(output_vcf_path, 'rt') as result_vcf:
                observed_lines = result_vcf.read().splitlines(keepends=False)

            observed_lines = [x for x in observed_lines if not x.startswith("#")]

            with open(self.swaps_expected_genotypes_vcf, 'rt') as expected_vcf:
                expected_lines = expected_vcf.read().splitlines(keepends=False)

            if observed_lines != expected_lines:
                for line in difflib.context_diff(expected_lines, observed_lines,
                                                 fromfile=self.swaps_expected_genotypes_vcf, tofile=output_vcf_path):
                    sys.stderr.write(line + "\n")

                with open("test_swaps_genotypes.vcf", "wt") as out_file:
                    for x in observed_lines:
                        print(x, file=out_file)

                raise Exception("Swaps VCF output doesn't match! If this is expected and new behavior, "
                                "cp test_swaps_genotypes.vcf %s" % self.swaps_expected_genotypes_vcf)

            output_json_path = os.path.join(output_dir, "variants.json.gz")
            with gzip.open(output_json_path, 'rt') as result_json:
                observed = json.load(result_json)

            with open(self.swaps_expected_variants_json, 'rt') as expected_json:
                expected = json.load(expected_json)

            # event ordering is not guaranteed
            observed = sorted(observed, key=lambda x: x["ID"])
            expected = sorted(expected, key=lambda x: x["ID"])

            # remove temp file information
            for o in observed:
                del o["graph"]["model_name"]
            for e in expected:
                if "model_name" in e["graph"]:
                    del e["graph"]["model_name"]

            observed_lines = json.dumps(observed, sort_keys=True, indent=4, separators=(',', ': '))
            expected_lines = json.dumps(expected, sort_keys=True, indent=4, separators=(',', ': '))
            if observed_lines != expected_lines:
                for line in difflib.context_diff(expected_lines.split("\n"), observed_lines.split("\n"),
                                                 fromfile=self.swaps_expected_variants_json, tofile=output_json_path):
                    sys.stderr.write(line + "\n")

                with open("test_swaps_variants.json", "wt") as out_file:
                    json.dump(observed, out_file, sort_keys=True, indent=4, separators=(',', ': '))

                raise Exception("Swaps test converted variants don't match! If this is expected and new behavior, "
                                "cp test_swaps_variants.json %s" % self.swaps_expected_variants_json)

    @unittest.skipIf(not GRMPY_INSTALL,
                     "No compiled/install path specified, please set the GRMPY_INSTALL variable to point to an installation.")
    @unittest.skipIf(not HG38, "No hg38 reference fasta file was found.")
    def test_multigrmpy_pg_het_ins(self):
        import multigrmpy

        with tempfile.TemporaryDirectory() as output_dir:
            args = ADict({
                "input": self.hg38_input_vcf,
                "manifest": self.hg38_manifest,
                "reference": self.hg38_reference,
                "output": output_dir,
                "grmpy": os.path.join(os.path.dirname(__file__), "module-wrapper.sh") + " " + os.path.join(
                    GRMPY_INSTALL, "bin", "grmpy"),
                "verbose": False,
                "quiet": True,
                "logfile": None,
                "write_alignments": False,
                "infer_read_haplotypes": False,
                "graph_sequence_matching": True,
                "klib_sequence_matching": False,
                "kmer_sequence_matching": False,
                "bad_align_uniq_kmer_len": 0,
                "threads": 1,
                "sample_threads": 1,
                "genotyping_parameters": None,
                "max_reads_per_event": 10000,
                "split_type": "lines",
                "read_length": 150,
                "max_ref_node_length": 1000,
                "ins_info_key": "SEQ",
                "graph_type": "alleles",
                "alt_splitting": True,
                "retrieve_reference_sequence": False,
                "scratch_dir": None,
                "keep_scratch": True,
            })
            output_json_path = os.path.join(output_dir, "genotypes.json.gz")
            multigrmpy.run(args)

            # compare genotyping results
            with gzip.open(output_json_path, 'rt') as result_json:
                observed = json.load(result_json)

            with open(self.hg38_expected_genotypes_json, 'rt') as expected_json:
                expected = json.load(expected_json)

            # uncomment to keep observed json
            # check if genotypes are the same
            observed_lines = json.dumps(observed, sort_keys=True, indent=4, separators=(',', ': '))
            expected_lines = json.dumps(expected, sort_keys=True, indent=4, separators=(',', ': '))
            if observed_lines != expected_lines:
                for line in difflib.context_diff(expected_lines.split("\n"), observed_lines.split("\n"),
                                                 fromfile=self.hg38_expected_genotypes_json, tofile=output_json_path):
                    sys.stderr.write(line + "\n")

                with open("test_gt.json", "wt") as out_file:
                    json.dump(observed, out_file, sort_keys=True, indent=4, separators=(',', ': '))

                raise Exception("Genotyping results don't match! If this is expected and new behavior, "
                                "cp test_gt.json %s" % self.hg38_expected_genotypes_json)

            # check if variants are the same
            output_json_path = os.path.join(output_dir, "variants.json.gz")
            with gzip.open(output_json_path, 'rt') as result_json:
                observed = json.load(result_json)

            # uncomment to keep observed json

            with open(self.hg38_expected_variants_json, 'rt') as expected_json:
                expected = json.load(expected_json)

            # contains temp file locations
            del observed[0]["graph"]["model_name"]
            if "model_name" in expected[0]["graph"]:
                del expected[0]["graph"]["model_name"]

            observed_lines = json.dumps(observed, sort_keys=True, indent=4, separators=(',', ': '))
            expected_lines = json.dumps(expected, sort_keys=True, indent=4, separators=(',', ': '))
            if observed_lines != expected_lines:
                for line in difflib.context_diff(expected_lines.split("\n"), observed_lines.split("\n"),
                                                 fromfile=self.hg38_expected_genotypes_json, tofile=output_json_path):
                    sys.stderr.write(line + "\n")

                with open("test_variants.json", "wt") as out_file:
                    json.dump(observed, out_file, sort_keys=True, indent=4, separators=(',', ': '))

                raise Exception("Converted variants don't match! If this is expected and new behavior, "
                                "cp test_variants.json %s" % self.hg38_expected_variants_json)


if __name__ == "__main__":
    unittest.main()
