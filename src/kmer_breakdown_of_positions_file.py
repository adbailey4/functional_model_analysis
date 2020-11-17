#!/usr/bin/env python
"""Write M kmer breakdown by position given a positions file and reference"""
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
from py3helpers.utils import time_it, all_string_permutations


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
    assert os.path.isdir(args.output_dir), "{} is not a directory".format(args.output_dir)
    assert os.path.exists(args.reference), "{} does not exist".format(args.reference)
    assert os.path.exists(args.positions_file), "{} does not exist".format(args.positions_file)

    positions_data = pd.read_csv(args.positions_file, names=["chr", "start", "strand", "find", "replace"], sep="\t")
    positions_data["kmer"] = np.nan
    # reference handler and reverse complement handler
    rh = ReferenceHandler(args.reference)
    rc = ReverseComplement()
    chromosome_data = {chromosome: rh.get_sequence(chromosome, 0, rh.get_chr_sequence_length(chromosome))
                       for chromosome in rh.fasta.references}

    alphabet = "ACGMT"
    kmer_length = 6

    def get_kmer(chromosome, pos, strand, replace=None):
        try:
            seq = chromosome_data[chromosome][(pos - kmer_length) + 1:pos + kmer_length]
            if strand == "-":
                seq = rc.reverse_complement(seq)
            if replace is not None:
                seq = seq[:kmer_length - 1] + replace + seq[kmer_length:]
        except Exception as e:
            print(e, chromosome, pos, strand)
        return seq

    mod_pos_data = positions_data.loc[positions_data['replace'] == "M"].copy()
    mod_pos_data.loc[:, "kmer"] = np.vectorize(get_kmer)(mod_pos_data['chr'],
                                                         mod_pos_data['start'],
                                                         mod_pos_data['strand'],
                                                         "M")

    kmers = {k: 0 for k in all_string_permutations(alphabet, kmer_length)}
    large_kmers = set(mod_pos_data['kmer'])
    for l_kmer in large_kmers:
        for i in range(kmer_length):
            k = l_kmer[i:kmer_length + i]
            if len(k) == kmer_length:
                kmers[k] += 1
    m_kmers = [x for x, y in kmers.items() if x.count("M") == 1]
    found_m_only_kmers = {x: y for x, y in kmers.items() if y > 0 and x.count("M") == 1}
    print(f"Number of M kmers: {len(m_kmers)}")
    print(f"Number of found M kmers: {len(found_m_only_kmers)}")

    c_pos_data = positions_data.loc[positions_data['replace'] == "C"].copy()
    c_pos_data.loc[:, 'kmer'] = np.vectorize(get_kmer)(c_pos_data['chr'],
                                                       c_pos_data['start'],
                                                       c_pos_data['strand'],
                                                       "C")
    filter_c_pos_data = c_pos_data[~c_pos_data["kmer"].str.contains('|'.join(["N", "W", "Y"]), regex=True)]

    kmers = {k: 0 for k in all_string_permutations(alphabet, kmer_length)}
    large_kmers = set(filter_c_pos_data['kmer'])
    for l_kmer in large_kmers:
        for i in range(kmer_length):
            k = l_kmer[i:kmer_length + i]
            if len(k) == kmer_length:
                kmers[k] += 1
    no_m_kmers = [x for x, y in kmers.items() if x.count("M") == 0 and x.count("C") > 0]
    found_no_m_kmers = {x: y for x, y in kmers.items() if y > 0 and x.count("M") == 0 and x.count("C") > 0}
    print(f"Number of Canonical kmers: {len(no_m_kmers)}")
    print(f"Number of found Canonical kmers: {len(found_no_m_kmers)}")


if __name__ == '__main__':
    print(time_it(main)[1])
