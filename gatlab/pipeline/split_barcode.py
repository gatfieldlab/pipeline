#!/usr/bin/env python3

"""
A simple barcode splitter for fastq formatted outputs
"""

import sys
import time
import itertools
import multiprocessing
import argparse
from collections import defaultdict
from xopen import xopen

__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2019, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.3.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"

VERSION = __version__


NTS = set(['A', 'T', 'C', 'G'])
PRIMERS = ("i7", "i5")
GZIP_THREADS = 12
SPLITTER_THREADS = 16
FILE_EXT = "_split.fastq.gz"

def get_args():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("barcode_file", metavar="barcode-file",
                        type=argparse.FileType('r'),
                        help="file containing sample names for barcode combinations."
                        "This file must have i) no header, ii) tab/space delimited "
                        "i7 i5 sample_name in each line")
    parser.add_argument("-g", "--gzip-threads", metavar="GZIP_T", type=int,
                        default=GZIP_THREADS, help="number of gzip (pigz) threads")
    parser.add_argument("-p", "--parser-threads", metavar="PARSER_T", type=int,
                        default=SPLITTER_THREADS, help="number of parser threads")
    parser.add_argument("-x", "--ext", default=FILE_EXT,
                        help="output file suffix/extension")
    parser.add_argument("-v", "--verbose", type=int, help="Verbosity level. "
                        "0 = silent, 1 = log, 2 = detailed log, 3+ = debug",
                        default=1)
    return parser.parse_args()

def get_barcodes(barcode_file):
    barcodes = {primer: set([]) for primer in PRIMERS}
    barcode_lens = {primer: None for primer in PRIMERS}
    sample_name = {}
    for line in barcode_file:
        barcode_seq, sname = line.strip().rsplit(maxsplit=1)
        barcode_seq = barcode_seq.split()
        assert len(barcode_seq) == len(PRIMERS), \
            "need exactly '{}' barcode sequences".format(len(PRIMERS))
        for primer, seq in zip(PRIMERS, barcode_seq):
            barcodes[primer].add(seq)
        sample_name[tuple(barcode_seq)] = sname
    barcode_len = {primer: set([len(seq) for seq in seqs])
                   for primer, seqs in barcodes.items()}
    for primer, seq_len in barcode_len.items():
        assert len(barcode_len[primer]) == 1, \
            "non unique sequence lengths for '{}'".format(primer)
        barcode_lens[primer] = seq_len.pop()
    return barcodes, barcode_lens, sample_name

def get_oneoff(barcode_seqs):
    one_off_seq = {}
    for seq in barcode_seqs:
        for pos_i, nt in enumerate(seq):
            for alt_nt in NTS - set([nt]):
                one_off_seq[seq[:pos_i] + alt_nt + seq[pos_i+1:]] = seq
    return one_off_seq

def mp_splitter(primers, barcodes, barcode_lens, sample_name,
                one_off_seq, verbose, job_q, result_q):
    i7_start_off = sum(barcode_lens.values()) + 2
    i7_end_off = barcode_lens['i5'] + 2
    i5_start_off = barcode_lens['i5'] + 1
    i5_end_off = 1
    for identifier, seq, sep, qual in iter(job_q.get, "STOP"):
        iden_len = len(identifier)
#        cur_seq = identifier.strip().rsplit(":", 1)[1].split("+", 1)
        cur_seq = [identifier[iden_len - i7_start_off:iden_len - i7_end_off],
                   identifier[iden_len - i5_start_off:iden_len - i5_end_off]]
        status = []
        cur_barcodes = []
        for i, primer in primers:
            if cur_seq[i] in barcodes[primer]:
                label = "perfect"
                cur_barcodes.append(cur_seq[i])
            elif cur_seq[i] in one_off_seq[primer]:
                if verbose > 1:
                    label = cur_seq[i]
                else:
                    label = "one_off"
                cur_barcodes.append(one_off_seq[primer][cur_seq[i]])
            else:
                label = "unident"
            status.append(primer + "_" + label)
        if len(cur_barcodes) == len(PRIMERS):
            sname = sample_name[tuple(cur_barcodes)]
            sbyte = (identifier+seq+sep+qual).encode()
        else:
            sname = None
            sbyte = None
        result_q.put((sname, sbyte, "__".join(status)))
    job_q.put("STOP") # for others
    result_q.put((None, None, "STOP")) # for manager

def mp_manager(writers, splitter_threads, verbose, result_q):
    counter = defaultdict(int)
    finished_splitters = 0
    while finished_splitters < splitter_threads:
        sname, sbyte, status = result_q.get()
        if status == "STOP":
            finished_splitters += 1
        else:
            if sname:
                writers[sname].put(sbyte)
            counter["Total sequences"] += 1
            counter[status] += 1
    for writer in writers.values():
        writer.put("STOP")
    if verbose:
        for sname, count in counter.items():
            sys.stderr.write("{}: {}\n".format(sname, count))

def mp_writer(sname, extension, gzip_threads, verbose, write_q):
    filename = "".join([sname, extension])
    out_handle = xopen(filename, "wb", threads=gzip_threads)
    check_point = time.perf_counter()
    wait_time = 0
    write_time = 0
    items = 0
    for sbyte in iter(write_q.get, "STOP"):
        wait_point = time.perf_counter()
        wait_time += wait_point - check_point
        out_handle.write(sbyte)
        check_point = time.perf_counter()
        write_time += check_point - wait_point
        items += 1
    if verbose >= 3:
        sys.stderr.write("This writer is finished: {}\n".format(sname))
        sys.stderr.write(
            "Mean wait: {} / mean write: {} / total written: {}\n".format(
                wait_time / items, write_time / items, items))
    out_handle.close()

def main():
    """
    Splits incoming fastq input based on barcode information (i7 and i5)
    1:N:0:ATCACG+ATAGAGGC
    """
    args = get_args()
    barcodes, barcode_lens, sample_name = get_barcodes(args.barcode_file)
    one_off_seq = {primer: get_oneoff(seqs) for primer, seqs in barcodes.items()}
    job_queue = multiprocessing.Queue()
    # Multiprocessing setup
    result_queue = multiprocessing.Queue()
    write_queues = {}
    splitters = []
    writers = []
    for sname in sample_name.values():
        writer_queue = multiprocessing.Queue()
        write_queues[sname] = writer_queue
        writer = multiprocessing.Process(
            target=mp_writer, args=(sname, args.ext, args.gzip_threads,
                                    args.verbose, writer_queue))
        writer.start()
        writers.append(writer)
    for i in range(args.parser_threads):
        splitter = multiprocessing.Process(
            target=mp_splitter,
            args=(list(enumerate(PRIMERS)), barcodes, barcode_lens, sample_name,
                  one_off_seq, args.verbose, job_queue, result_queue))
        splitter.start()
        splitters.append(splitter)
    manager = multiprocessing.Process(
        target=mp_manager, args=(write_queues, args.parser_threads,
                                 args.verbose, result_queue))
    manager.start()
    start_point = time.perf_counter_ns()
    if args.verbose:
        sys.stderr.write("split_barcode.py version {}\n".format(VERSION))
    for identifier, seq, sep, qual in itertools.zip_longest(*[sys.stdin]*4):
        job_queue.put((identifier, seq, sep, qual))
    job_queue.put("STOP")
    if args.verbose >= 3:
        sys.stderr.write("Already done with the job pipe. Time: {} sec\n".format(
            (time.perf_counter_ns() - start_point)/10**9))
    for splitter in splitters:
        splitter.join()
    if args.verbose >= 3:
        sys.stderr.write("All parser threads were joined. Time: {} sec\n".format(
            (time.perf_counter_ns() - start_point)/10**9))
    manager.join()
    if args.verbose >= 3:
        sys.stderr.write("Manager thread was joined. Time: {} sec\n".format(
            (time.perf_counter_ns() - start_point)/10**9))
    for writer in writers:
        writer.join()
    if args.verbose >= 3:
        sys.stderr.write("All writer threads were joined. Time: {} sec\n".format(
            (time.perf_counter_ns() - start_point)/10**9))
