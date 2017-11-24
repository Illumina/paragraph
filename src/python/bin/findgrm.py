#!/usr/bin/env python3

# coding=utf-8
#
# Copyright (c) 2010-2015 Illumina, Inc.
# All rights reserved.
#
# This file is distributed under the simplified BSD license.
# The full text can be found here (and in LICENSE.txt in the root folder of
# this distribution):
#
# https://github.com/Illumina/licenses/blob/master/Simplified-BSD-License.txt
#
# Script to find the installation path of grmpy s.t. things work both from a
# dev checkout and compiled / deployed version.
#

import os
import sys

_this_script_dir = os.path.abspath(os.path.dirname(__file__))

_dirs_to_try = [
    os.path.join(_this_script_dir, "..", "lib"),
    os.path.join(_this_script_dir, "..", "lib", "python3"),
]
_dirs_to_try = list(map(os.path.abspath, _dirs_to_try))

found = False
GRM_BASE = os.path.abspath(os.path.join(_this_script_dir, ".."))
for d in _dirs_to_try:
    if os.path.exists(os.path.join(d, "grm", "versionfallback.py")):
        sys.path.append(d)
        found = True
        break

if not found:
    raise Exception("Cannot find grmpy Python libraries in %s" % str(list(_dirs_to_try)))

try:
    import grm.version as version
except:  # pylint: disable=bare-except
    import grm.versionfallback as version


def main():
    print("Grmpy version      %s" % version.__version__)
    print("Grmpy base         %s" % GRM_BASE)
    print("Python sys.path    %s" % str(sys.path))


if __name__ == "__main__":
    main()
