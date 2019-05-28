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
import subprocess
import shutil

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


class TestParagraphPhasing(unittest.TestCase):
    def setUp(self):
        self.test_data_dir = os.path.join(
            GRMPY_ROOT, "share", "test-data", "paragraph", "phasing")

    @unittest.skipIf(not GRMPY_INSTALL, "No compiled/install path specified, please set the GRMPY_INSTALL variable to point to an installation.")
    def test_phasing(self):
        bam = os.path.join(self.test_data_dir, "na12878-many-haps.bam")
        in_json = os.path.join(self.test_data_dir, "long-phasing.json")
        output_file = tempfile.NamedTemporaryFile(delete=False, suffix=".json")
        output_file.close()

        subprocess.check_call([os.path.join(GRMPY_INSTALL, "bin", "paragraph"),
                               "-r", HG38,
                               "-b", bam,
                               "-g", in_json,
                               "-o", output_file.name,
                               "-E", "1",
                               "--bad-align-uniq-kmer-len", "64"])

        expected_fn = os.path.join(self.test_data_dir, "expected.json")
        with open(expected_fn, "rt") as fp:
            expected = json.load(fp)
            expected = [expected["phased_path_groups"],
                        expected["paths"]]

        with open(output_file.name, "rt") as fp:
            observed = json.load(fp)
            observed = [observed["phased_path_groups"],
                        observed["paths"]]

        if expected != observed:
            shutil.copy(output_file.name, "phasing_output.json")
            print("Difference between observed 'phasing_output.json' and expected ({})".format(expected_fn), file=sys.stderr)
        self.assertEqual(expected, observed)
        os.remove(output_file.name)


if __name__ == "__main__":
    unittest.main()
