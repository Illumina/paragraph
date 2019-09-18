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
import sys
import json
import unittest
import tempfile
import gzip

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


if __name__ == "__main__":
    unittest.main()
