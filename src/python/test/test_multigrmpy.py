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

GRMPY_ROOT = os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", ".."))

sys.path.append(os.path.join(GRMPY_ROOT,
                             "src", "python", "bin"))

GRMPY_INSTALL = None
if "GRMPY_INSTALL" in os.environ:
    GRMPY_INSTALL = os.environ["GRMPY_INSTALL"]


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

    @unittest.skipIf(not GRMPY_INSTALL, "No compiled/install path specified, please set the GRMPY_INSTALL variable to point to an installation.")
    def test_multigrmpy(self):
        import multigrmpy

        with tempfile.TemporaryDirectory() as output_dir:
            args = ADict({
                "input": self.input_vcf,
                "manifest": self.manifest,
                "reference": self.reference,
                "output": output_dir,
                "grmpy": os.path.join(os.path.dirname(__file__), "module-wrapper.sh") + " " + os.path.join(GRMPY_INSTALL, "bin", "grmpy"),
                "verbose": False,
                "quiet": True,
                "logfile": None,
                "threads": 1,
                "split_type": "lines",
                "read_length": 150,
                "max_ref_node_length": 1000,
                "graph_type": "alleles",
                "retrieve_reference_sequence": False,
                "scratch_dir": None,
                "keep_scratch": True,
                "use_em": False
            })
            output_json_path = os.path.join(output_dir, "genotypes.json.gz")
            multigrmpy.run(args)
            with gzip.open(output_json_path, 'rt') as result_json:
                observed = json.load(result_json)
                self.assertEqual(observed[0]["samples"]["sample1"]["gt"]["GT"], "S1/S1")
                self.assertEqual(observed[1]["samples"]["sample2"]["gt"]["GT"], "S1/S1")


if __name__ == "__main__":
    unittest.main()
