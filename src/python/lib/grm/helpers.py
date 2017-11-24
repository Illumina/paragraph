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
# August 2017
#
# Helper functions
#
# Author:
#
# Peter Krusche <pkrusche@illumina.com>
#

import logging


def parse_region(region_string):
    """
    Parse a region string into chrom, start, end
    :param region_string: string like "chrom:start-end"
    :type region_string: str
    :return: chrom, start, end
    """
    region_string = region_string.replace(",", "")
    sp = region_string.split(":", 1)
    chrom = sp[0]
    start = None
    end = None
    if len(sp) > 1:
        startend = sp[1].split("-", 1)
        start = int(startend[0])
        if len(startend) > 1:
            end = int(startend[1])
    return chrom, start, end


class LoggingWriter(object):  # pylint: disable=too-few-public-methods
    """ Helper class to write tracebacks to log file
    """

    def __init__(self, level):
        self.level = level

    def write(self, message):
        message = message.replace("\n", "")
        if message:
            logging.log(self.level, message)
