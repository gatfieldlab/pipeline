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
__version__ = "0.2.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"

VERSION = __version__
DEF_CUTOFF = 10000


def main():
    """
    Merges gzipped fastq files reading one read from each at a time
    """
    parser = argparse.ArgumentParser()
    parser.add_argument(files, help="gzipped fastq files to process")
    parser.add_argument("--max-seq", help="max number of sequences to merge",
                        type=int, default=DEF_CUTOFF)
    args = parser.parse_args()
    handlers = [gzip.open(seqfile, 'rt') for seqfile in args.files]
    iterators = [itertools.zip_longest(*[handle]*4) for handle in handlers]
    print(handlers, file=sys.stderr)
    count = 0
    while iterators:
        ditch_list = []
        for i, iterator in enumerate(iterators):
            try:
                print(*next(iterator), sep="", end="")
                count += 1
            except StopIteration:
                ditch_list.append(i)
                print("appended to ditchlist:", ditch_list)
        for di in ditch_list[::-1]:
            iterators.pop(di)
            handlers.pop(di).close()
        if count >= args.max_seq:
            break
    print("Total =", count, file=sys.stderr)

if __name__ == "__main__":
    sys.exit(main())
