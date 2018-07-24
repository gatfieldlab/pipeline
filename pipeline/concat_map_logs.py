#!/usr/bin/env python3

"""
Doc string goes to here
"""

import os
import re
import sys


READS_REGEX = re.compile(r'^(\d+) reads')
ALIGN_REGEX = re.compile(r'(\d+) .+ aligned (.+)$')


def get_bw_logfiles(log_path, base_name, extension="bw2.log"):
    """
    Doc goes here
    """
    file_regex = re.compile('^{}\.(.+)\.{}'.format(base_name, extension))
    file_list = os.listdir(log_path)
    log_files = []
    for file in file_list:
        file_match = file_regex.match(file)
        if file_match:
            log_type = file_match.group(1)
            log_files.append((log_type, os.path.join(log_path, file)))
    return log_files


def get_star_log(log_path, base_name):
    """
    Doc goes here
    """
    log_file = os.path.join(log_path, base_name, "Log.final.out")
    with open(log_file) as open_file:
        for line in open_file:
            try:
                key, value = tuple([val.strip() for val
                                    in line.strip().split("|")])
            except ValueError:
                continue
            if key == "Uniquely mapped reads number":
                uniq_map = int(value)
            elif key == "Number of reads mapped to multiple loci":
                multi_map = int(value)
            elif key == "Number of reads mapped to too many loci":
                super_map = int(value)
            elif key == "Number of input reads":
                all_reads = int(value)
        try:
            total_mapped = uniq_map + multi_map + super_map
            unmapped = all_reads - total_mapped
        except NameError as err:
            print("Error while parsing {}".format(log_file))
            print(err)
            sys.exit(1)
        return (all_reads, total_mapped, unmapped)


def main():
    """
    Doc goes here
    """
    log_path = sys.argv[2]
    star_log_path = sys.argv[3]
    with open(sys.argv[1]) as samples:
        for sample in samples:
            sample = sample.strip()
            star_log = get_star_log(star_log_path, sample)
            log_files = get_bw_logfiles(log_path, sample)
            for log_type, log_file in log_files:
                align_stat = {}
                with open(log_file) as open_file:
                    for line in open_file:
                        read_match = READS_REGEX.match(line)
                        if read_match:
                            tot_reads = read_match.group(1)
                        align_match = ALIGN_REGEX.search(line)
                        if align_match:
                            align_stat[align_match.group(2)] = align_match.group(1)
                    align_stat['total_aligned'] = str(int(tot_reads) - int(align_stat['0 times']))
                    keys = ['total_aligned']
                    for t in keys:
                        print("\t".join([sample, log_type, t, align_stat[t]]))
                    if log_type == 'mouse_cDNA':
                        assert int(align_stat['0 times']) == star_log[0]
            print("\t".join([sample, "genomic", "total_aligned", str(star_log[1])]))
            print("\t".join([sample, "unmapped", "total_aligned", str(star_log[2])]))
