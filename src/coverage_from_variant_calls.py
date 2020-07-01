#!/usr/bin/env python
"""Get kmer coverage from variant call file and reference"""
########################################################################
# File: coverage_from_variant_calls.py
#  executable: coverage_from_variant_calls.py
#
# Author: Andrew Bailey
# History: 06/30/20 Created
########################################################################

import os
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from signalalign.hiddenMarkovModel import HmmModel
from py3helpers.utils import list_dir, all_string_permutations, time_it
from py3helpers.seq_tools import ReferenceHandler, ReverseComplement


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--variants', '-v', action='store', nargs='+',
                        dest='variants', required=True, type=str, default=None,
                        help="Path to variant call file")
    parser.add_argument('--output_dir', '-o', action='store',
                        dest='output_dir', required=True, type=str, default=None,
                        help="Path to directory where kmer count data will be stored")
    parser.add_argument('--reference', '-r', action='store',
                        dest='reference', required=True, type=str, default=None,
                        help="Path to reference sequence")
    parser.add_argument('--rna', action='store_true',
                        dest='rna', required=False,
                        help="Boolean option if models are for rna")
    parser.add_argument('--alphabet', action='store',
                        dest='alphabet', required=True, type=str,
                        help="alphabet")
    parser.add_argument('--kmer_length', action='store',
                        dest='kmer_length', required=True, type=int,
                        help="kmer_len")

    args = parser.parse_args()
    return args


def main():
    args = parse_args()
    print(args)
    assert os.path.exists(args.reference), "{} does not exist".format(args.reference)
    assert os.path.isdir(args.output_dir), "{} is not a directory".format(args.output_dir)
    for variant_path in args.variants:
        assert os.path.exists(variant_path), "{} does not exist".format(variant_path)
    rh = ReferenceHandler(args.reference)
    rc = ReverseComplement()
    kmers = {k: 0 for k in all_string_permutations(args.alphabet, args.kmer_length)}
    for variant_path in args.variants:
        print(variant_path)
        local_kmers = {k: 0 for k in all_string_permutations(args.alphabet, args.kmer_length)}

        variant_data = pd.read_csv(variant_path)
    #         read_id            028a34d4-2a7a-44e7-ab23-305915996ec8
    #         contig                                          RDN18-1
    #         reference_index                                     973
    #         strand                                                +
    #         variants                                             Aa
    #         prob1                                          0.986967
    #         prob2                                          0.013033
    #         prob3                                               NaN

        for i, r in variant_data.iterrows():
            seq = rh.get_sequence(r.contig, (r.reference_index-args.kmer_length)+1, r.reference_index+args.kmer_length)
            if r.strand == "-":
                seq = rc.reverse_complement(seq)
            for i in range(args.kmer_length):
                kmers[seq[i:args.kmer_length+i]] += 1
                local_kmers[seq[i:args.kmer_length+i]] += 1

        total_zeros = 0
        for x, y in local_kmers.items():
            if y == 0:
                total_zeros += 1
        print("Kmers Covered: {}/{}".format(len(local_kmers)-total_zeros, len(local_kmers)))
        print("Average coverage: {}".format(np.average(list(local_kmers.values()))))
    total_zeros = 0
    for x, y in kmers.items():
        if y == 0:
            total_zeros += 1

    print("TOTAL Kmers Covered: {}/{}".format(len(kmers)-total_zeros, len(kmers)))
    print("TOTAL Average coverage: {}".format(np.average(list(kmers.values()))))


if __name__ == '__main__':
    print(time_it(main)[1])
