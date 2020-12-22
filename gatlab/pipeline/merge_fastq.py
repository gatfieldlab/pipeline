#!/usr/bin/env python3

"""
A simple pipe that combines a number of fastq files by sequentially extracting
single reads from each file
"""

import sys
import itertools
import gzip
import argparse

__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2020, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.3.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"

VERSION = __version__
DEF_CUTOFF = 0


def main():
    """
    Merges gzipped fastq files reading one sequence from each at a time
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("files", help="gzipped fastq files to process",
                        nargs="+")
    parser.add_argument(
        "--max-seq", help="max number of sequences to merge, 0 implies no limit",
        type=int, default=DEF_CUTOFF)
    args = parser.parse_args()
    handlers = [gzip.open(seqfile, 'rt') for seqfile in args.files]
    iterators = [itertools.zip_longest(*[handle]*4) for handle in handlers]
    print("[merge_fastq]", [handle.name for handle in handlers],
          "are opened to be merged.", file=sys.stderr)
    count = 0
    while iterators:
        ditch_list = []
        for i, iterator in enumerate(iterators):
            try:
                print(*next(iterator), sep="", end="")
                count += 1
            except StopIteration:
                ditch_list.append(i)
            except BrokenPipeError:
                print("[merge_fastq] Broken pipe detected - processing stopped",
                      file=sys.stderr)
                break
        for di in ditch_list[::-1]:
            iterators.pop(di)
            handle = handlers.pop(di)
            print("[merge_fastq]", handle.name, "finished, closing",
                  file=sys.stderr)
            handle.close()
        if args.max_seq and count >= args.max_seq:
            break
    for handle in handlers:
        handle.close()
    print("[merge_fastq]", args.files, "Total =", count, file=sys.stderr)
