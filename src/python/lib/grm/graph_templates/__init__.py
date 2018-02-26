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
# October 2017
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import os
from . import shortdeletion
from . import longdeletion
from . import shortswap
from . import longswap
from . import insertion


def make_graph(reference, event):
    """
    Create graph dictionary from event description
    :param reference: filename for reference fasta
    :param event: event description
                  {
                    "chrom": ...
                    "start": ...
                    "end": ...

                    // optional
                    "ins": ...

                    // optional
                    "flank": ...
                  }
    :return: type, graph dictionary that can be written to JSON for processing in paragraph
    """
    assert os.path.exists(reference)
    try:
        flank = event["flank"]
    except KeyError:
        flank = 150

    try:
        ins = event["ins"]
    except KeyError:
        ins = ""

    # needs to be deletion, insertion or both
    ref_len = event["end"] - event["start"] + 1
    is_del = ref_len > 0
    assert is_del or len(ins) > 0  # pylint: disable=len-as-condition

    chrom = event["chrom"]
    start = min(event["start"], event["end"])
    end = max(event["start"], event["end"])

    if is_del and len(ins) == 0:   # pylint: disable=len-as-condition
        if ref_len <= 2 * flank:
            return "del", shortdeletion.make_graph(chrom, start, end, flank)
        else:
            return "longdel", longdeletion.make_graph(chrom, start, end, flank)
    elif is_del:  # -> len(ins) > 0
        if ref_len <= 2 * flank:
            return "swap", shortswap.make_graph(chrom, start, end, ins, flank)
        else:
            return "longswap", longswap.make_graph(chrom, start, end, ins, flank)
    else:  # not is_del and len(ins) > 0
        return "ins", insertion.make_graph(chrom, start, ins, flank)
