#!/usr/bin/env python
"""Get top N positions covered by reads"""
########################################################################
# File: top_n_most_covered_positions.py
#  executable: top_n_most_covered_positions.py
#
# Author: Andrew Bailey
# History: 04/20/21 Created
########################################################################

import os
import sys
import pandas as pd
import numpy as np
from argparse import ArgumentParser
import pysam

from py3helpers.utils import time_it, list_dir


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--bam', '-b', action='store',
                        dest='bam', required=True, type=str, default=None,
                        help="Path to bam file")
    parser.add_argument('--positions_file', '-p', action='store',
                        dest='positions_file', required=False, type=str, default=None,
                        help="Path to positions_file")
    parser.add_argument('--output_dir', '-o', action='store',
                        dest='output_dir', required=True, type=str, default=None,
                        help="Path to output directory")
    parser.add_argument('-n', action='store',
                        dest='n', required=True, type=int, default=100,
                        help="number of positions to write")
    parser.add_argument('--positions_directory', '-d', action='store',
                        dest='positions_directory', required=False, type=str, default=None,
                        help="Path to positions_directory")

    args = parser.parse_args()
    return args


class WriteTopN(object):
    def __init__(self, bam_file):
        self.bam_file = bam_file
        self.alignment_handle = pysam.AlignmentFile(self.bam_file, "rb")

    def get_coverage(self, chr, start):
        return np.sum([x[0] for x in self.alignment_handle.count_coverage(chr, start, start+1, quality_threshold=0)])

    def write_top_n_covered_positions_file(self, positions_file, output_dir, n):
        output_path = os.path.join(output_dir, f"top_{n}_covered_{os.path.basename(positions_file)}")
        cov_out_dir = os.path.join(output_dir, "coverage")
        if not os.path.exists(cov_out_dir):
            os.mkdir(cov_out_dir)
        coverage_out_path = os.path.join(cov_out_dir, f"top_{n}_covered_{os.path.basename(positions_file)}.coverage")

        positions_data = pd.read_csv(positions_file, names=["chr", "start", "strand", "find", "replace"], sep="\t")

        coverages = np.vectorize(self.get_coverage)(positions_data["chr"],
                                               positions_data['start'])
        positions_data["coverage"] = coverages
        positions_data = positions_data[positions_data["coverage"] < 40]
        sorted_pos = positions_data.sort_values("coverage", ascending=False)
        write_out = sorted_pos[:n]

        write_out.to_csv(output_path, sep="\t", index=False,
                         columns=["chr", "start", "strand", "find", "replace"], header=False)
        write_out.to_csv(coverage_out_path, sep="\t", index=False,
                         columns=["chr", "start", "strand", "find", "replace", "coverage"], header=False)


def main():
    args = parse_args()
    assert os.path.isdir(args.output_dir), "{} is not a directory".format(args.output_dir)
    assert os.path.exists(args.bam), "{} does not exist".format(args.bam)

    if args.positions_directory is not None:
        assert os.path.exists(args.positions_directory), "{} does not exist".format(args.positions_file)
    else:
        assert args.positions_file is not None, "Must pass in --positions_file or --positions_directory"
        assert os.path.exists(args.positions_file), "{} does not exist".format(args.positions_file)
    output_dir = args.output_dir
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    bam = args.bam
    top_n = WriteTopN(bam_file=bam)
    n = args.n

    if args.positions_directory is not None:
        for positions_file in list_dir(args.positions_directory, ext="positions"):
            print(f"Running on {positions_file}")
            top_n.write_top_n_covered_positions_file(positions_file, output_dir, n)
    else:
        print(f"Running on {args.positions_file}")
        top_n.write_top_n_covered_positions_file(args.positions_file, output_dir, n)


if __name__ == '__main__':
    print(time_it(main)[1])
