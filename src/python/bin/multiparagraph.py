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
import argparse
import subprocess
import tempfile
import multiprocessing
import json
import logging
import gzip
import traceback
import itertools
import pipes

import findgrm  # pylint: disable=unused-import

from grm.helpers import LoggingWriter
from grm.graph_templates import make_graph


def run_paragraph(event_and_args):
    """ Wrapper function to run paragraph
    :param event_and_args: tuple of event to validate and args structure
    :return:
    """
    event = event_and_args[0]
    args = event_and_args[1]
    tempfiles = []
    exception = ""
    error_log = ""
    error = False
    try:
        # paragraph input
        tf = tempfile.NamedTemporaryFile(dir=args.scratch_dir, mode="wt", suffix=".json", delete=False)
        tempfiles.append(tf.name)

        if "graph" not in event:
            event["type"], event["graph"] = make_graph(args.ref, event)
            logging.info("%s:%i type=%s len=%i id=%i", event["chrom"],
                         event["start"], event["type"], event["len"], event["n_ev"])
        elif "type" not in event:
            event["type"] = "custom"

        json.dump(event["graph"], tf)
        tf.close()

        # paragraph output
        output_name = tf.name + ".output.json"
        tempfiles.append(output_name)
        error_log = tf.name + ".output.log"
        tempfiles.append(error_log)

        commandline = args.paragraph
        commandline += " -r %s" % pipes.quote(args.ref)
        commandline += " -b %s" % pipes.quote(args.bam)
        commandline += " -g %s" % pipes.quote(tf.name)
        commandline += " -o %s" % pipes.quote(output_name)
        commandline += " --log-file %s" % pipes.quote(error_log)
        commandline += " --threads %i" % args.paragraph_threads

        if args.extended_output:
            commandline += " -E 1"

        event["commandline"] = commandline
        o = subprocess.check_output(commandline, shell=True, stderr=subprocess.STDOUT)

        try:
            o = o.decode("utf-8")
        except:  # pylint: disable=bare-except
            o = str(o)

        for line in o.split("\n"):
            if line:
                logging.warning(line)

        with open(output_name, "rt") as f:
            event["graph"] = json.load(f)

    except Exception:  # pylint: disable=broad-except
        logging.error("Exception when running paragraph on %s", str(event))
        logging.error('-' * 60)
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        logging.error('-' * 60)
        exception = traceback.format_exc()
        error = True
    except BaseException:  # pylint: disable=broad-except
        logging.error("Exception when running paragraph on %s", str(event))
        logging.error('-' * 60)
        traceback.print_exc(file=LoggingWriter(logging.ERROR))
        logging.error('-' * 60)
        exception = traceback.format_exc()
        error = True
    finally:
        if error:
            try:
                with open(error_log, "rt") as f:
                    for line in f:
                        line = str(line).strip()
                        logging.error(line)
            except:  # pylint: disable=bare-except
                pass
            event["error"] = {
                "exception": exception
            }
        if not args.keep_scratch:
            for f in tempfiles:
                try:
                    os.remove(f)
                except:  # pylint: disable=bare-except
                    pass
    return event


def make_argument_parser():
    """
    :return: an argument parser
    """
    parser = argparse.ArgumentParser("Multiparagraph")

    parser.add_argument("input", nargs="+")

    parser.add_argument("-b", "--bam", help="BAM file name",
                        type=str, dest="bam", required=True)

    parser.add_argument("-o", "--output", help="Output file name",
                        type=str, dest="output", required=True)

    parser.add_argument("-r", "--reference-sequence", help="Reference FASTA",
                        type=str, dest="ref", required=True)

    parser.add_argument("-E", "--extended-output", help="Run paragraph with -E 1",
                        dest="extended_output", action="store_true", default=False)

    parser.add_argument("--max-events", dest="max_events", type=int, default=None,
                        help="Only do the first n events.")

    parser.add_argument("--min-length", dest="min_length", type=int, default=0,
                        help="Minimum event length.")

    parser.add_argument("--event-threads", "-t", dest="threads", type=int, default=multiprocessing.cpu_count(),
                        help="Number of events to process in parallel.")

    parser.add_argument("--paragraph-threads", "-T", dest="paragraph_threads", type=int, default=1,
                        help="Number of threads for parallel read processing.")

    parser.add_argument("--keep-scratch", dest="keep_scratch", default=None, action="store_true",
                        help="Do not delete temp files.")

    parser.add_argument("--scratch-dir", dest="scratch_dir", default=None,
                        help="Directory for temp files")

    parser.add_argument("--paragraph", dest="paragraph", default=os.path.join(findgrm.GRM_BASE, "bin", "paragraph"),
                        help="Path to the paragraph executable")

    parser.add_argument("--logfile", dest="logfile", default=None,
                        help="Write logging information into file rather than to stderr")

    verbosity_options = parser.add_mutually_exclusive_group(required=False)

    verbosity_options.add_argument("--verbose", dest="verbose", default=False, action="store_true",
                                   help="Raise logging level from warning to info.")

    verbosity_options.add_argument("--quiet", dest="quiet", default=False, action="store_true",
                                   help="Set logging level to output errors only.")

    return parser


def run(args):
    """ Run wrapper for testing """
    if args.verbose:
        loglevel = logging.INFO
    elif args.quiet:
        loglevel = logging.ERROR
    else:
        loglevel = logging.WARNING

        # reinitialize logging
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    logging.basicConfig(filename=args.logfile,
                        format='%(asctime)s %(levelname)-8s %(message)s',
                        level=loglevel)

    assert os.path.exists(args.bam)

    all_events = []
    for filename in args.input:
        if filename.endswith(".gz"):
            f = gzip.open(filename, "rt")
        else:
            f = open(filename, "rt")

        events = json.load(f)
        f.close()

        if not isinstance(events, list):
            raise Exception("Input JSON must contain a list of events.")

        for e in events:
            if not isinstance(e, dict):
                raise Exception("Invalid event description: %s" % str(e))
            if "ins" in e:
                e["ins_len"] = len(e["ins"])
            else:
                e["ins_len"] = 0

            e["del_len"] = max(0, e["end"] - e["start"] + 1)
            e["len"] = max(e["ins_len"], e["del_len"])
            if "samples" in e:
                del e["samples"]

        all_events += events

    if args.max_events is not None:
        all_events = all_events[:args.max_events]

    all_events = [e for e in all_events if e["len"] >= args.min_length]
    i = 0
    for e in all_events:
        e["n_ev"] = i
        i += 1

    logging.info("Number of events: %i", len(all_events))

    with multiprocessing.Pool(args.threads) as pool:
        results = pool.map(run_paragraph, zip(all_events, itertools.repeat(args)))

    if args.output.endswith(".gz"):
        of = gzip.open(args.output, "wt")
        json.dump(results, of)
    else:
        of = open(args.output, "wt")
        json.dump(results, of, sort_keys=True, indent=4, separators=(',', ': '))

    of.close()


def main():
    parser = make_argument_parser()
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
