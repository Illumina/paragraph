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


class TestMultiParagraph(unittest.TestCase):
    def setUp(self):
        self.test_data_dir = os.path.join(
            GRMPY_ROOT, "share", "test-data", "multiparagraph")

    @unittest.skipIf(not GRMPY_INSTALL, "No compiled/install path specified, please set the GRMPY_INSTALL variable to point to an installation.")
    def test_multiparagraph(self):
        bam = os.path.join(self.test_data_dir, "reads.bam")
        ref = os.path.join(self.test_data_dir, "dummy.fa")
        candidates = os.path.join(self.test_data_dir, "candidates.json")
        output_file = tempfile.NamedTemporaryFile(delete=False, suffix=".json")
        output_file.close()

        args = ADict({
            "bam": bam,
            "ref": ref,
            "input": [candidates],
            "output": output_file.name,
            "verbose": False,
            "quiet": True,
            "logfile": None,
            "max_events": None,
            "min_length": 0,
            "extended_output": False,
            "paragraph_threads": 1,
            "threads": 1,
            "scratch_dir": None,
            "keep_scratch": False,
            "paragraph": os.path.join(os.path.dirname(__file__), "module-wrapper.sh") + " " + os.path.join(GRMPY_INSTALL, "bin", "paragraph")
        })
        import multiparagraph
        multiparagraph.run(args)

        with open(os.path.join(self.test_data_dir, "expected.json"), "rt") as f:
            expected = json.load(f)

        with open(output_file.name, "rt") as f:
            observed = json.load(f)

        for x in expected:
            del x["graph"]["bam"]
            del x["graph"]["reference"]
            del x["commandline"]
        for x in observed:
            del x["graph"]["bam"]
            del x["graph"]["reference"]
            del x["commandline"]

        expected = json.dumps(expected,
                              sort_keys=True, indent=4, separators=(',', ': ')).splitlines()
        observed = json.dumps(observed,
                              sort_keys=True, indent=4, separators=(',', ': ')).splitlines()

        if expected != observed:
            d = difflib.Differ()
            diff = d.compare(expected, observed)
            print("\n".join(diff))
        self.assertEqual(expected, observed)
        os.remove(output_file.name)


if __name__ == "__main__":
    unittest.main()
