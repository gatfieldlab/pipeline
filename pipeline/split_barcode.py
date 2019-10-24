#!/usr/bin/env python3

"""
A simple barcode splitter for fastq formatted outputs
"""

import sys
import itertools
from collections import defaultdict

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
    Splits incoming fastq input based on barcode informatione (i5 and i7)
    1:N:0:ATCACG+ATAGAGGC
    """

    total_count = 0
    split_count = defaultdict(int)
    output_base = "sample"
    sample_name = {}
    one_off_seq = {}
    i5s = set([])
    i7s = set([])
    i5_range = list(range(8))
    nts = set(['A', 'T', 'C', 'G'])
    with open("barcodes.txt") as f:
        for line in f:
            i7, i5, sname = line.strip().split()
            i5s.add(i5)
            i7s.add(i7)
            sample_name[(i7, i5)] = sname
    for i5 in i5s:
        for pos_i, nt in enumerate(i5):
            for alt_nt in nts - set([nt]):
                one_off_seq[i5[:pos_i] + alt_nt + i5[pos_i+1:]] = i5
    openfiles = {sname: open("_".join([output_base, sname, "split.fastq"]), "w") for
                             sname in sample_name.values()}
    # def _get_one_off(seq):
    #     res = []
    #     for i5 in i5s:
    #         if sum([i5[i] == seq[i] for i in i5_range]) == 1:
    #             res.append(i5)
    #     return res
    for identifier, seq, sep, qual in itertools.zip_longest(*[sys.stdin]*4):
        total_count += 1
        i7, i5 = identifier.strip().rsplit(":", 1)[1].split("+", 1)
        split_ok = False
        if (i7, i5) in sample_name:
            sname = sample_name[(i7, i5)]
            cname = sname
            split_ok = True
            split_count["perfect"] += 1
        else:
            # one_off_candidates = _get_one_off(i5)
            if i7 in i7s and i5 in one_off_seq:
                corr_i5 = one_off_seq[i5]
                sname = sample_name[(i7, corr_i5)]
                cname = "_".join(["oneoff", i7, corr_i5])
                split_count["one_off"] += 1
            else:
                sname = None
                cname = "_".join([i7, i5])
                split_count["unidentified"] += 1
        if split_ok:
            openfiles[sname].write(identifier+seq+sep+qual)            
        split_count[cname] += 1
    for outfile in openfiles.values():
        outfile.close()
    sys.stderr.write("split_barcode.py version {}\n".format(VERSION))
    sys.stderr.write("Total number of sequences: {}\n".format(total_count))
    for sname, count in split_count.items():
        sys.stderr.write("{}: {}\n".format(sname, count))
        
#    sys.stderr.write("Number of sequences filtered out: {}\n".format(filtered_count))
#    sys.stderr.write("Filter used: [min:{}, max:{}]\n".format(min_seq_len, max_seq_len))
