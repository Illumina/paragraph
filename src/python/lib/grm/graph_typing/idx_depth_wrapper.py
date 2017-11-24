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
# September 2017
#
# idx depth related function for runGraphTyping.py
#
# Author:
#
# Sai Chen <schen6@illumina.com>
#

import subprocess
import pipes
import os
import logging

import findgrm


def generate_idxdepth_json(bam_manifest,
                           reference,
                           idx_depth_dir,
                           threads=1,
                           idxdepth_binary=os.path.join(findgrm.GRM_BASE, "bin", "idxdepth")):
    """
    Generate idxdepth JSON for all samples

    :param bam_manifest: list of BAM file information records
    :type bam_manifest: list of grm.graph_typing.sample_manifest.ManifestElement
    :param reference: reference fasta file name
    :param idx_depth_dir: output directory.
    :param threads: number of threads
    :param idxdepth_binary: custom binary for indexdepth
    """
    for bam_element in bam_manifest:
        idx_name = bam_element.name + ".idxdepth"
        idx_log = idx_name + ".log"
        outpath = os.path.join(idx_depth_dir, idx_name)
        logpath = os.path.join(
            idx_depth_dir, idx_log)
        try:
            idx_cmd = idxdepth_binary
            idx_cmd += " -r %s" % pipes.quote(reference)
            idx_cmd += " -b %s" % pipes.quote(bam_element.path)
            idx_cmd += " -o %s" % pipes.quote(outpath)
            idx_cmd += " --threads %i" % threads
            idx_cmd += " --log-file %s" % pipes.quote(logpath)
            subprocess.check_call(idx_cmd, shell=True)
        except:  # pylint: disable=bare-except
            logging.error(
                "Idxdepth failed for %s, see %s for details", outpath, logpath)
            raise
