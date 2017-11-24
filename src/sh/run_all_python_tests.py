#! /usr/bin/env python3

import os
import sys
import unittest


def main():
    project_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "python"))

    testsuite = unittest.TestLoader().discover(os.path.join(project_dir, "test"))
    lib_results = unittest.TextTestRunner(verbosity=3).run(testsuite)

    print("%d tests run in test/, %d failures, %d errors" %
          (lib_results.testsRun, len(lib_results.failures), len(lib_results.errors)))

    if len(lib_results.failures) > 0:
        sys.exit(1)

    if len(lib_results.errors) > 0:
        sys.exit(1)


if __name__ == '__main__':
    main()
