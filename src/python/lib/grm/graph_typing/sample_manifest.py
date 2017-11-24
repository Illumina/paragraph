#!/usr/bin/env python3
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
#
# October 2017
#
# functions related to manifest data structure
#
#
# Author:
#
# Sai Chen <schen6@illumina.com>
#
#
import errno
import os
import json


class ManifestElement(object):
    """
    record and modification of one manifest record
    """

    def __init__(self):
        self.name = ""
        self.path = ""
        self.depth = 0
        self.read_len = 0

    def load(self, line):
        fields = line.split('\t')
        if len(fields) != 2 and len(fields) != 4:
            raise Exception(
                "Manifest must have either 2 columns or 4 columns.")
        self.name = fields[0]
        self.path = fields[1]
        if len(fields) == 4:
            if "." in fields[3]:
                raise Exception("Read length is not an integer.")
            try:
                self.depth = fields[2]
                self.read_len = fields[3]
            except ValueError:
                raise Exception(
                    "Depth and read length fields should only contain numbers.")

    def append_to_path(self, extra_path):
        self.path = os.path.join(extra_path, self.path)

    def missing_stats(self):
        return bool(self.depth == 0 or self.read_len == 0)

    def to_string(self):
        out_str = self.name + "\t" + self.path + "\t" + \
            str(self.depth) + "\t" + str(self.read_len)
        return out_str


def load_manifest(manifest_path, use_relative_manifest_path=False):
    """
    load manifest from file. One line for one record
    """
    if use_relative_manifest_path:
        extra_path = os.path.dirname(manifest_path)
    manifest = []
    sample_names = set()
    with open(manifest_path, "rt") as manifest_file:
        for line in manifest_file:
            if line.startswith('#'):
                continue
            line = line.rstrip()
            element = ManifestElement()
            element.load(line)

            if use_relative_manifest_path:
                element.append_to_path(extra_path)
            elif not os.path.exists(element.path):
                # check paths relative to manifest
                element.append_to_path(os.path.dirname(
                    os.path.abspath(manifest_path)))

            # check if file exists
            if not os.path.exists(element.path):
                raise FileNotFoundError(
                    errno.ENOENT, os.strerror(errno.ENOENT), element.path)

            # duplicate sample names are bad because we name folders after them
            if element.name in sample_names:
                raise Exception(
                    "Duplicate sample name in manifest: %s" % element.name)
            sample_names.add(element.name)

            manifest.append(element)
    return manifest


def generate_and_update_manifest(idx_depth_dir, manifest, new_manifest_path):
    """
    update manifest data from generated idxdepth results
    dump manifest into the new path
    """
    for element in manifest:
        sample_idx_path = os.path.join(
            idx_depth_dir, element.name + ".idxdepth")
        if not os.path.isfile(sample_idx_path):
            message = "Missing idxdepth result at sample: " + sample_idx_path
            raise Exception(message)
        with open(sample_idx_path, 'r') as f:
            js = json.load(f)
            try:
                depth = js["autosome"]["depth"]
                read_len = js["read_length"]
            except KeyError:
                message = "Missing required fields in idxdepth output at sample: " + element.sample
                raise Exception(message)
            element.depth = depth
            element.read_len = read_len
            f.close()
    with open(new_manifest_path, 'w') as f:
        f.write("#ID\tPath\tDepth\tRead_len\n")
        for element in manifest:
            f.write(element.to_string() + "\n")
        f.close()


def count_missing_stats(manifest):
    """
    count the number of samples in manifest with missing depth or read length
    """
    num_missing = 0
    for element in manifest:
        if element.missing_stats():
            num_missing += 1
    return num_missing
