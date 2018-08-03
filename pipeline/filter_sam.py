#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Filters SAM output
"""

import sys
from argparse import ArgumentParser
import fileinput
from signal import signal, SIGPIPE, SIG_DFL
from accessories import utils


__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2013-2018, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "2.1.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"


signal(SIGPIPE, SIG_DFL)



def get_args():
    parser = ArgumentParser()
    parser.add_argument(
        'files', nargs='*', help='specify input SAM files')
    parser.add_argument(
        '-o', '--output', help='specify output file. The default is stdout')
    parser.add_argument(
        '-c', '--command', help='', type=str, default='filter',
        choices=['filter', 'passtru', 'basic', 'rrna'])
    parser.add_argument(
        '-v', '--verbose', help='a little more verbosity', action='store_true')
    parser.add_argument(
        '-s', '--size-file', help='size (bp) file for references in SAM')
    parser.add_argument(
        '-k', '--keep-filtered', action='store_true',
        help='keep the filtered alignments in a separate SAM file')
    parser.add_argument(
        '-f', '--file-kept', type=str, default='filtered.sam',
        help='filename for --keep-filtered. if empty default is "filtered.sam"')
    parser.add_argument(
        '--version', help='prints version and quits', action='store_true')
    return parser.parse_args()


def process_cur_reads(cur_reads, line, command, keep_filtered, kept_file, rrna_sizes):
    if command == 'passtru':
        sys.stdout.write(line)
        return True
    parsed = line.strip().split('\t')
    bit_op = int(parsed[1])
    if bit_op & 0x4:
        if keep_filtered:
            kept_file.write(line)
        return cur_reads
    if command == 'basic':
        sys.stdout.write(line)
        return True
    # command == 'filter' or 'rrna'
    read_reverse = bool(bit_op & 0x10)
    if read_reverse and command == 'filter':
        if keep_filtered:
            kept_file.write(line)
        return cur_reads
    if command == 'rrna':
        try:
            ref_name = parsed[2]
            ref_len = rrna_sizes[ref_name]
        except KeyError as err:
            sys.stderr.write(
                "SAM seems to include a references for which the size "
                "information is not available.\n{}".format(err.message))
            sys.exit(1)
    else:
        ref_name = None
        ref_len = None
    read_id = parsed[0]
    op_as = None
    for f in parsed[11:]:
        op, t, val = f.split(':')
        if op == 'AS':
            op_as = int(val)
            break
    assert op_as != None, "NO 'AS' tag for this read???"
    if (cur_reads and read_id != cur_reads['read_id']) or not cur_reads:
        if cur_reads:
            if command == 'filter':
                sys.stdout.write(cur_reads['lines'])
            if command == 'rrna':
                hist = cur_reads['hist']
                if len(hist) > 1:
                    size_priority = sorted(range(len(hist)), key=lambda k: hist[k], reverse=True)
                    max_len = hist[size_priority[0]][0]
                    curlines = cur_reads['lines'].split('\n')
                    # Size filtered only -- no directions at this point
                    for refindex in size_priority:
                        if hist[refindex][0] == max_len:
                            sys.stdout.write("{}\n".format(curlines[refindex]))
                        else:
                            if keep_filtered:
                                kept_file.write("{}\n".format(curlines[refindex]))
                            else:
                                break
                else:
                    sys.stdout.write(cur_reads['lines'])
        cur_reads = {'read_id': read_id, 'AS': op_as, 'hist': [(ref_len, ref_name, read_reverse)],
                     'lines': line} # first hit is always kept if not 0x4 (x10)
    else: # cur_reads set and read_id is same -- anything else ??
        if op_as == cur_reads['AS']:
            cur_reads['lines'] += line
            if command == 'rrna':
                cur_reads['hist'].append((ref_len, ref_name, read_reverse))
        else:
            if keep_filtered:
                kept_file.write(line)
    return cur_reads


def main():
    args = get_args()
    if args.version:
        sys.stdout.write('filter_sam version {}\n'.format(__version__))
        sys.exit(0)

    if args.output and args.output != '-':
        sys.stdout = open(args.output, 'w')
    if args.verbose:
        sys.stderr.write('Opened {} for output\n'.format(sys.stdout.name))
    if args.command == 'passtru':
        args.keep_filtered = False
    if args.keep_filtered:
        kept_file = open(args.file_kept, 'w')
        if args.verbose:
            sys.stderr.write('Opened {} for keeping filtered alignments\n'.format(args.file_kept))
    rrna_sizes = {}
    if args.command == 'rrna':
        if not args.size_file:
            sys.stderr.write("Can't filter rRNA without the reference size information.")
            sys.exit(1)
        else:
            with open(args.size_file) as sizef:
                for line in sizef:
                    parsed = line.strip().split('\t')
                    rrna_sizes[parsed[0]] = int(parsed[1])
            if args.verbose:
                sys.stderr.write('Extracted length for {} references from {}\n'.format(len(rrna_sizes), args.size_file))

    with utils.measureTime('Finished filtering {}'.format(','.join(args.files))):
        samfiles = fileinput.input(args.files)
        cur_reads = None
        while True:
            try:
                line = next(samfiles)
            except StopIteration:
                if args.command == 'filter' and cur_reads:
                    sys.stdout.write(cur_reads['lines'])
                break
            if line[0] == '@':
                if not cur_reads:
                    sys.stdout.write(line)
                if fileinput.isfirstline():
                    sys.stderr.write("Processing input file: '{}'\n".format(fileinput.filename()))
                continue
            else:
                cur_reads = process_cur_reads(cur_reads, line, args.command,
                                              args.keep_filtered, kept_file, rrna_sizes)
