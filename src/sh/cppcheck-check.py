#!/usr/bin/env python
# Simple checker for whether cppcheck found errors

from __future__ import print_function
import sys
import xml.etree.ElementTree as ElementTree

e = ElementTree.parse(sys.argv[1])


errors = e.findall('errors')

errorcount = 0
for es in errors:
    for ei in es.findall('error'):
        print("Error: cppcheck:%s" % ei.attrib["id"])
        errorcount += 1


sys.exit(0 if not errorcount else 1)
