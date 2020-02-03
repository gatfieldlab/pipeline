#!/usr/bin/env python3

"""
A simple size filter for fastq formatted outputs
"""

import sys
import itertools


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2015-2019, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.2.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"

VERSION = __version__

def main():
    """
    Filters incoming fastq input based read size (MIN, MAX)
    """
    if len(sys.argv) > 1 and sys.argv[1] == '-v':
        sys.stdout.write("consume {}\n".format(VERSION))
        sys.exit(0)
    try:
        min_seq_len = int(sys.argv[1])
    except IndexError:
        min_seq_len = 21 #for RB, 26
    try:
        max_seq_len = int(sys.argv[2])
    except IndexError:
        max_seq_len = 60 #for RB, 35

    total_count = 0
    filtered_count = 0

    # to account for newline character
    eff_min_seq_len = min_seq_len + 1
    eff_max_seq_len = max_seq_len + 1
    if len(sys.argv) > 3:
        outfile = open(sys.argv[3], 'w')
    else:
        outfile = sys.stdout

    for identifier, seq, sep, qual in itertools.zip_longest(*[sys.stdin]*4):
        total_count += 1
        if len(seq) < eff_min_seq_len or len(seq) > eff_max_seq_len:
            filtered_count += 1
        else:
            outfile.write(identifier+seq+sep+qual)
    outfile.close()
    sys.stderr.write("size_filter.py version {}\n".format(VERSION))
    sys.stderr.write("Total number of sequences: {}\n".format(total_count))
    sys.stderr.write("Number of sequences filtered out: {}\n".format(filtered_count))
    sys.stderr.write("Filter used: [min:{}, max:{}]\n".format(min_seq_len, max_seq_len))
