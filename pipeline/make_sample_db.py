#!/usr/bin/env python3

"""
utility to extract a tab-delimited flat text file to
serve as a simple database of samples for pipeline,
from files in an NGS output folder
"""

__author__ = "Bulak Arpat"
__copyright__ = "Copyright 2017-2019, Bulak Arpat"
__license__ = "GPLv3"
__version__ = "0.2.0"
__maintainer__ = "Bulak Arpat"
__email__ = "Bulak.Arpat@unil.ch"
__status__ = "Development"

import sys

def main():
    for line in sys.stdin:
        filename = line.strip()
        base, ext = filename.split('.', 1)
        try:
            sample, run_id, lane_id, file_id = base.rsplit('_', 3)
        except ValueError:
            sys.stderr.write('Failed parsing: {}\n'.format(base))
            continue
        annot = sample.split('_')
        sample_run_id = '_'.join([sample, run_id, lane_id])
        sys.stdout.write('\t'.join([run_id, lane_id, file_id,
                                    sample_run_id, base, filename] + annot) + '\n')
