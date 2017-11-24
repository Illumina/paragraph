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
# Test the graph SV workflow using the round trip examples
#

import os
import unittest
import tempfile
import shutil
import subprocess
import pipes

GRMPY_ROOT = os.path.abspath(os.path.join(
    os.path.dirname(__file__), ".."))

GRMPY_INSTALL = None
if "GRMPY_INSTALL" in os.environ:
    GRMPY_INSTALL = os.environ["GRMPY_INSTALL"]


class TestMain(unittest.TestCase):
    def setUp(self):
        self.out_name = tempfile.mkdtemp()
        hg19_locations = [
            "/Users/pkrusche/workspace/human_genome/hg19.fa",
            "/illumina/sync/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa",
            "/Users/schen6/Documents/hg19/genome.fa"
        ]
        if "HG19" in os.environ and os.environ["HG19"]:
            hg19_locations.insert(0, os.environ["HG19"])
        self.hg19 = None
        for f in hg19_locations:
            if os.path.exists(f) and os.path.isfile(f):
                self.hg19 = f
                break

        if not self.hg19:
            raise Exception(
                "hg19 fasta file not found in expected locations: %s" % str(hg19_locations))

        self.run_graph_typing = "python3 " + \
                                str(GRMPY_INSTALL) + "/bin/runGraphTyping.py"
        self.grm_install = None

    def tearDown(self):
        shutil.rmtree(self.out_name)

    @unittest.skipIf(not GRMPY_INSTALL, "No compiled/install path specified, please set the GRMPY_INSTALL variable to point to an installation.")
    def test_round_trip_json_input(self):
        input_json = GRMPY_INSTALL + "/share/test-data/round-trip-genotyping/candidates.json"
        manifest_path = GRMPY_INSTALL + "/share/test-data/round-trip-genotyping/samples.txt"
        ref_path = GRMPY_INSTALL + "/share/test-data/round-trip-genotyping/dummy.fa"
        out_path = tempfile.mkdtemp(suffix="test_json")
        cmd = self.run_graph_typing
        cmd += " -i %s" % pipes.quote(input_json)
        cmd += " -r %s" % pipes.quote(ref_path)
        cmd += " -m %s" % pipes.quote(manifest_path)
        cmd += " -o %s" % pipes.quote(out_path)
        cmd += " --relative-manifest-path"
        cmd += " --grmpy-binary " + pipes.quote(
            os.path.join(os.path.dirname(__file__), "module-wrapper.sh") + " " + os.path.join(GRMPY_INSTALL, "bin", "grmpy"))
        cmd += " --paragraph-binary " + pipes.quote(
            os.path.join(os.path.dirname(__file__), "module-wrapper.sh") + " " + os.path.join(GRMPY_INSTALL, "bin", "paragraph"))
        cmd += " --idxdepth-binary " + pipes.quote(
            os.path.join(os.path.dirname(__file__), "module-wrapper.sh") + " " + os.path.join(GRMPY_INSTALL, "bin", "idxdepth"))
        print(cmd)
        subprocess.check_call(cmd, shell=True)
        shutil.rmtree(out_path)

    @unittest.skipIf(not GRMPY_INSTALL, "No compiled/install path specified, please set the GRMPY_INSTALL variable to point to an installation.")
    def test_round_trip_vcf_input(self):
        input_vcf = GRMPY_INSTALL + "/share/test-data/round-trip-genotyping/candidates.vcf"
        manifest_path = GRMPY_INSTALL + "/share/test-data/round-trip-genotyping/samples.txt"
        ref_path = GRMPY_INSTALL + "/share/test-data/round-trip-genotyping/dummy.fa"
        out_path = tempfile.mkdtemp(suffix="test_json")
        cmd = self.run_graph_typing
        cmd += " -i %s" % pipes.quote(input_vcf)
        cmd += " -r %s" % pipes.quote(ref_path)
        cmd += " -m %s" % pipes.quote(manifest_path)
        cmd += " -o %s" % pipes.quote(out_path)
        cmd += " --relative-manifest-path"
        cmd += " --grmpy-binary " + pipes.quote(
            os.path.join(os.path.dirname(__file__), "module-wrapper.sh") + " " + os.path.join(GRMPY_INSTALL, "bin", "grmpy"))
        cmd += " --paragraph-binary " + pipes.quote(
            os.path.join(os.path.dirname(__file__), "module-wrapper.sh") + " " + os.path.join(GRMPY_INSTALL, "bin", "paragraph"))
        cmd += " --idxdepth-binary " + pipes.quote(
            os.path.join(os.path.dirname(__file__), "module-wrapper.sh") + " " + os.path.join(GRMPY_INSTALL, "bin", "idxdepth"))
        print(cmd)
        subprocess.check_call(cmd, shell=True)
        shutil.rmtree(out_path)

    @unittest.skipIf(not GRMPY_INSTALL, "No compiled/install path specified, please set the GRMPY_INSTALL variable to point to an installation.")
    def test_genotypes_vcf_input(self):
        input_vcf = GRMPY_INSTALL + "/share/test-data/genotyping_test_2/swaps.vcf"
        manifest_path = GRMPY_INSTALL + "/share/test-data/genotyping_test_2/samples.txt"
        ref_path = GRMPY_INSTALL + "/share/test-data/genotyping_test_2/swaps.fa"
        out_path = tempfile.mkdtemp(suffix="test_json")
        cmd = self.run_graph_typing
        cmd += " -i %s" % pipes.quote(input_vcf)
        cmd += " -r %s" % pipes.quote(ref_path)
        cmd += " -m %s" % pipes.quote(manifest_path)
        cmd += " -o %s" % pipes.quote(out_path)
        cmd += " --grmpy-binary " + pipes.quote(
            os.path.join(os.path.dirname(__file__), "module-wrapper.sh") + " " + os.path.join(GRMPY_INSTALL, "bin", "grmpy"))
        cmd += " --paragraph-binary " + pipes.quote(
            os.path.join(os.path.dirname(__file__), "module-wrapper.sh") + " " + os.path.join(GRMPY_INSTALL, "bin", "paragraph"))
        cmd += " --idxdepth-binary " + pipes.quote(
            os.path.join(os.path.dirname(__file__), "module-wrapper.sh") + " " + os.path.join(GRMPY_INSTALL, "bin", "idxdepth"))
        cmd += " --relative-manifest-path"
        print(cmd)
        subprocess.check_call(cmd, shell=True)

        subprocess.check_call("gunzip %s" % pipes.quote(os.path.join(out_path, "SWAPS.genotype.json.gz")),
                              shell=True)

        subprocess.check_call("diff -I '\"ID\": \"Allele graph' %s %s" % (
            pipes.quote(os.path.join(out_path, "SWAPS.genotype.json")),
            pipes.quote(GRMPY_INSTALL + "/share/test-data/genotyping_test_2/expected-genotypes.json")),
                              shell=True)

        shutil.rmtree(out_path)


if __name__ == '__main__':
    unittest.main()
