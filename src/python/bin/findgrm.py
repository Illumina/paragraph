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
