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
import sys
import os
import unittest
import subprocess
import tempfile
from pipes import quote

sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/../lib")
from grm.vcfgraph import graphUtils, variants, graphContainer  # pylint: disable=C0413

GRMPY_ROOT = os.path.abspath(os.path.join(
    os.path.dirname(__file__), "..", "..", ".."))


class AddVariantsToplevelTests(unittest.TestCase):
    def setUp(self):
        self.test_data_dir = os.path.join(GRMPY_ROOT, "share", "test-data", "paragraph")
        self.test_graphs = [self.test_data_dir + f"/variants/{f}.json" for f in ["ref"]]
        self.addVariants = os.path.join(GRMPY_ROOT, "src", "python", "bin", "addVariants.py")

    def test_examples(self):
        for x in self.test_graphs:
            print("Testing output for {}".format(x))
            tf1 = tempfile.NamedTemporaryFile(suffix=".json")
            tf1.close()
            try:
                expected_json = x.replace(".json", "-vars.json")
                subprocess.check_call(
                    f"python3 {quote(self.addVariants)} {quote(x)} {tf1.name}",
                    shell=True)
                subprocess.check_call(
                    f"python3 -mjson.tool --sort-keys {tf1.name} > {tf1.name + '.pp.json'}",
                    shell=True)
                subprocess.check_call(
                    f"diff {tf1.name + '.pp.json'} {expected_json}",
                    stdout=subprocess.DEVNULL, shell=True)
            except:
                from shutil import copy
                current_dir = os.path.abspath(os.path.dirname(__file__))
                copy(tf1.name + ".pp.json", os.path.join(current_dir, "test-failed.json"))
                os.chmod(os.path.join(current_dir, "test-failed.json"), 0o777)
                print("Failed output saved in",
                      os.path.join(current_dir, "test-failed.json"))
                raise
            finally:
                os.remove(tf1.name)
                os.remove(tf1.name + ".pp.json")


class AddVariantsUnitTests(unittest.TestCase):
    def test_ref_snv(self):
        graph = graphContainer.GraphContainer()
        n = graph.add_refNode("chr", 10, 20)
        var = {"start": 2, "end": 2, "alt": "C"}
        variants.add_variants(graph, {n['name']: [var]})
        nodeNames = [n['name'] for n in graph.nodes.values()]
        self.assertCountEqual(nodeNames, ["ref-chr:10-11", "ref-chr:12-12", "ref-chr:13-20", "chr:12-12:C"])
        self.assertEqual(graph.nodes["chr:12-12:C"]["sequence"], "C")
        left = graph.nodes["ref-chr:10-11"]
        right = graph.nodes["ref-chr:13-20"]
        alt = graph.nodes["chr:12-12:C"]
        ref = graph.nodes["ref-chr:12-12"]
        self.assertTrue(graph.has_edge(left, ref))
        self.assertTrue(graph.has_edge(left, alt))
        self.assertTrue(graph.has_edge(ref, right))
        self.assertTrue(graph.has_edge(alt, right))
        self.assertFalse(graph.has_edge(left, right))
        self.assertFalse(graph.has_edge(ref, alt))

    def test_alt_snv(self):
        graph = graphContainer.GraphContainer()
        n = graph.add_altNode("chr", 10, 20, "ATCGATCG")
        var = {"start": 2, "end": 2, "alt": "T"}
        variants.add_variants(graph, {n['name']: [var]})
        nodeNames = [n['name'] for n in graph.nodes.values()]
        self.assertCountEqual(nodeNames, ["chr:10-11:AT", "chr:12-12:C", "chr:13-20:GATCG", "chr:12-12:T"])
        self.assertEqual(graph.nodes["chr:10-11:AT"]["sequence"], "AT")
        self.assertEqual(graph.nodes["chr:13-20:GATCG"]["sequence"], "GATCG")
        self.assertEqual(graph.nodes["chr:12-12:C"]["sequence"], "C")
        self.assertEqual(graph.nodes["chr:12-12:T"]["sequence"], "T")

    def test_insertion(self):
        graph = graphContainer.GraphContainer()
        n = graph.add_altNode("chr", 10, 17, "ATCGATCG")
        var = {"start": 3, "end": 2, "alt": "TTT"}
        variants.add_variants(graph, {n['name']: [var]})
        nodeNames = [n['name'] for n in graph.nodes.values()]
        self.assertCountEqual(nodeNames, ["chr:10-12:ATC", "chr:13-17:GATCG", "chr:13-12:TTT"])
        left = graph.nodes["chr:10-12:ATC"]
        right = graph.nodes["chr:13-17:GATCG"]
        ins = graph.nodes["chr:13-12:TTT"]
        self.assertTrue(graph.has_edge(left, right))
        self.assertTrue(graph.has_edge(left, ins))
        self.assertTrue(graph.has_edge(ins, right))
        self.assertFalse(graph.has_edge(ins, left))
        self.assertFalse(graph.has_edge(right, ins))

    def test_deletion(self):
        graph = graphContainer.GraphContainer()
        n = graph.add_altNode("chr", 10, 17, "ATCGATCG")
        var = {"start": 2, "end": 4, "alt": ""}
        variants.add_variants(graph, {n['name']: [var]})
        graphUtils.remove_empty_nodes(graph)
        nodeNames = [n['name'] for n in graph.nodes.values()]
        self.assertCountEqual(nodeNames, ["chr:10-11:AT", "chr:12-14:CGA", "chr:15-17:TCG"])
        left = graph.nodes["chr:10-11:AT"]
        right = graph.nodes["chr:15-17:TCG"]
        ins = graph.nodes["chr:12-14:CGA"]
        self.assertTrue(graph.has_edge(left, right))
        self.assertTrue(graph.has_edge(left, ins))
        self.assertTrue(graph.has_edge(ins, right))
        self.assertFalse(graph.has_edge(ins, left))
        self.assertFalse(graph.has_edge(right, ins))

    def test_var_begin(self):
        graph = graphContainer.GraphContainer()
        r = graph.add_refNode("chr", 1, 9)
        n = graph.add_altNode("chr", 10, 17, "ATCGATCG")
        graph.add_edge(r, n, ["foo"])
        var = {"start": 0, "end": 0, "alt": "G"}
        variants.add_variants(graph, {n['name']: [var]})
        graphUtils.remove_empty_nodes(graph)
        left = graph.nodes["ref-chr:1-9"]
        right = graph.nodes["chr:11-17:TCGATCG"]
        ref = graph.nodes["chr:10-10:A"]
        alt = graph.nodes["chr:10-10:G"]
        self.assertEqual(len(graph.nodes), 4)
        self.assertTrue(graph.has_edge(left, ref))
        self.assertTrue(graph.has_edge(left, alt))
        self.assertTrue(graph.has_edge(ref, right))
        self.assertTrue(graph.has_edge(alt, right))
        self.assertFalse(graph.has_edge(left, right))
        self.assertFalse(graph.has_edge(ref, alt))
        self.assertCountEqual(graph.get_edge(left['name'], ref['name'])['sequences'], ["foo"])
        self.assertCountEqual(graph.get_edge(left['name'], alt['name'])['sequences'], ["foo"])

    def test_var_end(self):
        graph = graphContainer.GraphContainer()
        r = graph.add_refNode("chr", 18, 20)
        n = graph.add_altNode("chr", 10, 17, "ATCGATCG")
        graph.add_edge(n, r, ["foo"])
        var = {"start": 7, "end": 7, "alt": "C"}
        variants.add_variants(graph, {n['name']: [var]})
        graphUtils.remove_empty_nodes(graph)
        left = graph.nodes["chr:10-16:ATCGATC"]
        right = graph.nodes["ref-chr:18-20"]
        ref = graph.nodes["chr:17-17:G"]
        alt = graph.nodes["chr:17-17:C"]
        self.assertEqual(len(graph.nodes), 4)
        self.assertTrue(graph.has_edge(left, ref))
        self.assertTrue(graph.has_edge(left, alt))
        self.assertTrue(graph.has_edge(ref, right))
        self.assertTrue(graph.has_edge(alt, right))
        self.assertFalse(graph.has_edge(left, right))
        self.assertFalse(graph.has_edge(ref, alt))
        self.assertCountEqual(graph.get_edge(ref['name'], right['name'])['sequences'], ["foo"])
        self.assertCountEqual(graph.get_edge(alt['name'], right['name'])['sequences'], ["foo"])

    def test_ins_end(self):
        graph = graphContainer.GraphContainer()
        r = graph.add_refNode("chr", 18, 20)
        n = graph.add_altNode("chr", 10, 17, "ATCGATCG")
        graph.add_edge(n, r, ["foo"])
        var = {"start": 8, "end": 7, "alt": "CCC"}
        variants.add_variants(graph, {n['name']: [var]})
        graphUtils.remove_empty_nodes(graph)
        left = graph.nodes["chr:10-17:ATCGATCG"]
        right = graph.nodes["ref-chr:18-20"]
        alt = graph.nodes["chr:18-17:CCC"]
        self.assertEqual(len(graph.nodes), 3)
        self.assertTrue(graph.has_edge(left, right))
        self.assertTrue(graph.has_edge(left, alt))
        self.assertTrue(graph.has_edge(alt, right))
        self.assertFalse(graph.has_edge(alt, left))
        self.assertFalse(graph.has_edge(right, alt))
        self.assertCountEqual(graph.get_edge(left['name'], right['name'])['sequences'], ["foo"])
        self.assertCountEqual(graph.get_edge(alt['name'], right['name'])['sequences'], ["foo"])

    def test_overlapping_deletion(self):
        graph = graphContainer.GraphContainer()
        n = graph.add_altNode("chr", 10, 17, "ATCGATCG")
        varDel = {"start": 2, "end": 4, "alt": ""}
        varSNV = {"start": 4, "end": 4, "alt": "C"}
        variants.add_variants(graph, {n['name']: [varDel, varSNV]})
        graphUtils.remove_empty_nodes(graph)
        nodeNames = [n['name'] for n in graph.nodes.values()]
        self.assertCountEqual(nodeNames,
                              ["chr:10-11:AT", "chr:12-13:CG", "chr:14-14:A", "chr:14-14:C", "chr:15-17:TCG"])
        left = graph.nodes["chr:10-11:AT"]
        right = graph.nodes["chr:15-17:TCG"]
        alt = graph.nodes["chr:14-14:C"]
        ref = graph.nodes["chr:14-14:A"]
        ins = graph.nodes["chr:12-13:CG"]
        self.assertTrue(graph.has_edge(left, right))
        self.assertTrue(graph.has_edge(left, ins))
        self.assertTrue(graph.has_edge(ins, ref))
        self.assertTrue(graph.has_edge(ref, right))
        self.assertTrue(graph.has_edge(ins, alt))
        self.assertTrue(graph.has_edge(alt, right))
        self.assertFalse(graph.has_edge(left, ref))
        self.assertFalse(graph.has_edge(left, alt))
        self.assertFalse(graph.has_edge(ref, ins))


if __name__ == "__main__":
    unittest.main()
