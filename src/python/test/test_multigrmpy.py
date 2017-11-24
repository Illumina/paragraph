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
import difflib
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
        self.test_pg_data_dir = os.path.join(
            GRMPY_ROOT, "share", "test-data", "multiparagraph")
        self.test_grmpy_data_dir = os.path.join(
            GRMPY_ROOT, "share", "test-data", "multigrmpy")

        self.tempfiles = []

    def tearDown(self):
        for x in self.tempfiles:
            try:
                os.remove(x)
            except:  # pylint: disable=bare-except
                pass

    @unittest.skipIf(not GRMPY_INSTALL, "No compiled/install path specified, please set the GRMPY_INSTALL variable to point to an installation.")
    def test_multigrmpy(self):
        ref = os.path.join(self.test_pg_data_dir, "dummy.fa")
        with tempfile.TemporaryDirectory() as output_dir:
            fo = tempfile.NamedTemporaryFile(
                mode="wt", delete=False, suffix=".manifest")
            self.tempfiles.append(fo.name)
            fo.write("#ID\tPath\tDepth\tRead_len\n")
            fo.write("Dummy\t" + os.path.join(self.test_pg_data_dir,
                                              "expected.json") + "\t1\t50\n")
            fo.close()

            args = ADict({
                "sample_name": "Dummy",
                "manifest": fo.name,
                "reference": ref,
                "output": output_dir,
                "id_identifier": "desc",
                "grmpy": os.path.join(os.path.dirname(__file__), "module-wrapper.sh") + " " + os.path.join(GRMPY_INSTALL, "bin", "grmpy"),
                "genotype_error_rate": 0.01,
                "min_overlap_bases": 16,
                "max_read_times": 10,
                "verbose": False,
                "quiet": True,
                "logfile": None,
                "threads": 1,
                "scratch_dir": None,
                "keep_scratch": False,
                "use_em": False
            })
            import multigrmpy
            multigrmpy.run(args)

            with open(os.path.join(self.test_grmpy_data_dir, "expected.json"), "rt") as f:
                expected = json.load(f)

            output_path = os.path.join(output_dir, "Dummy.genotype.json.gz")
            with gzip.open(output_path, "rt") as f:
                observed = json.load(f)

            expected = json.dumps(expected,
                                  sort_keys=True, indent=4, separators=(',', ': ')).splitlines()
            observed = json.dumps(observed,
                                  sort_keys=True, indent=4, separators=(',', ': ')).splitlines()

            if expected != observed:
                d = difflib.Differ()
                diff = d.compare(expected, observed)
                print("\n".join(diff))

            self.assertEqual(expected, observed)


if __name__ == "__main__":
    unittest.main()
