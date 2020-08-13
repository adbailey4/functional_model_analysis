#!/usr/bin/env python
"""Write kmer breakdown by position given a positions file and reference"""
########################################################################
# File: kmer_breakdown_of_positions_file.py
#  executable: kmer_breakdown_of_positions_file.py
#
# Author: Andrew Bailey
# History: 05/24/20 Created
########################################################################

import os
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from py3helpers.seq_tools import ReferenceHandler, ReverseComplement


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--reference', '-r', action='store',
                        dest='reference', required=True, type=str, default=None,
                        help="Path to reference")
    parser.add_argument('--positions_file', '-p', action='store',
                        dest='positions_file', required=True, type=str, default=None,
                        help="Path to positions_file")
    parser.add_argument('--output_dir', '-o', action='store',
                        dest='output_dir', required=True, type=str, default=None,
                        help="Path to output directory")

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    print(args)
    assert os.path.isdir(args.output_dir), "{} is not a directory".format(args.output_dir)
    assert os.path.exists(args.reference), "{} does not exist".format(args.reference)
    assert os.path.exists(args.positions_file), "{} does not exist".format(args.positions_file)

    positions_data = pd.read_csv(args.positions_file, names=["chr", "start", "strand", "find", "replace"], sep="\t")

    # reference handler and reverse complement handler
    rh = ReferenceHandler(args.reference)
    rc = ReverseComplement()
    chromosome_data = {chromosome: rh.get_sequence(chromosome, 0, rh.get_chr_sequence_length(chromosome))
                       for chromosome in rh.fasta.references}

    def get_adjacent_base(chromosome, pos, strand):
        try:
            base = chromosome_data[chromosome][pos+1]
            if strand == "-":
                base = chromosome_data[chromosome][pos-1]
                return rc.complement(base)
        except Exception as e:
            print(e, chromosome, pos, strand)
        return base

    def get_kmer(chromosome, pos, strand):
        try:
            base = chromosome_data[chromosome][pos-5:pos+6]
            if strand == "-":
                base = chromosome_data[chromosome][pos-5:pos+6]
                return rc.complement(base)
        except Exception as e:
            print(e, chromosome, pos, strand)
        return base

    positions_data['next_base'] = np.vectorize(get_kmer)(positions_data['chr'], positions_data['start'], positions_data['strand'])
    # cx_positions_data['kmer_base'] = np.vectorize(get_kmer_base)(cx_positions_data['chr'], cx_positions_data['start'], cx_positions_data['strand'])
    #
    # canonical_cx_positions_data['next_base'] = np.vectorize(get_adjacent_base)(canonical_cx_positions_data['chr'], canonical_cx_positions_data['start'], canonical_cx_positions_data['strand'])
    # cx_positions_data[cx_positions_data['next_base'] == "T"].to_csv("/home/ubuntu/bisulfite_methylation_analysis/positions/CT_methyl.positions", sep="\t", index=False, header=False, columns=["chr", "start", "strand", "find", "replace"])
    positions_data.to_csv(os.path.join(args.output_dir, os.path.basename(os.path.splitext(args.positions_file)[0])),
                          sep="\t", index=False, header=False,
                          columns=["chr", "start", "strand", "find", "replace", "next_base"])
