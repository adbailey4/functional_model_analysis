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
import time

import pandas as pd
import numpy as np
from py3helpers.utils import list_dir, all_string_permutations, time_it
from py3helpers.seq_tools import ReferenceHandler, ReverseComplement


def main():
    OUTPUT_DIR = "/home/ubuntu/mount/kmer_counts"
    positions_data = False
    keys = ["contig", "reference_index", "strand"]
    # RNA canonical
    REFERENCE = "/home/ubuntu/mount/download/RNA_rel2/reference/gencode.v27.transcripts.fa"
    VARIANT_HOME_DIRS = ["/home/ubuntu/mount/OICR_runs/all_runs/", "/home/ubuntu/mount/UBC_runs/all_runs/"]
    # VARIANT_HOME_DIRS = ["/home/ubuntu/mount/OICR_runs/test/", "/home/ubuntu/mount/OICR_runs/test/"]
    VARIANT_NAMES = ["/variant_calls/na12878_OICR_RNA_canonical.csv", "/variant_calls/na12878_UBC_RNA_canonical.csv"]
    ALPHABET = "ATGC"
    KMER_LENGTH = 5
    NAMES = ["OICR", "UBC"]

    # DNA canonical
    REFERENCE = "/home/ubuntu/bisulfite_methylation_analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    VARIANT_HOME_DIRS = ["/home/ubuntu/mount/FAB39088_runs/canonical_calling/all_runs/",
                         "/home/ubuntu/mount/FAF01169_runs/canonical_calling/all_runs/"]
    VARIANT_NAMES = ["/variant_calls/variant_calls.csv", "/variant_calls/variant_calls.csv"]
    ALPHABET = "ATGC"
    KMER_LENGTH = 6
    NAMES = ["FAB39088_canonical", "FAF01169_canonical"]

    # DNA mod
    # REFERENCE = "/home/ubuntu/bisulfite_methylation_analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
    # VARIANT_HOME_DIRS = ["/home/ubuntu/mount/FAB39088_runs/cpg_calling/all_runs/",
    #                      "/home/ubuntu/mount/FAF01169_runs/cpg_calling/all_runs/"]
    # NAMES = ["FAB39088_methyl", "FAF01169_methyl"]
    # VARIANT_NAMES = ["/variant_calls/variant_calls.csv", "/variant_calls/variant_calls.csv"]
    # ALPHABET = "ATGCM"
    # KMER_LENGTH = 6
    # POSITIONS_FILE = "/home/ubuntu/bisulfite_methylation_analysis/positions/all_mC.positions"
    # positions_data = pd.read_csv(POSITIONS_FILE, names=["contig", "reference_index", "strand", "find", "replace"],
    #                              sep="\t")
    if positions_data is not False:
        i2 = positions_data.set_index(keys).index

    assert os.path.exists(REFERENCE), "{} does not exist".format(REFERENCE)
    assert os.path.isdir(OUTPUT_DIR), "{} is not a directory".format(OUTPUT_DIR)
    rh = ReferenceHandler(REFERENCE)
    rc = ReverseComplement()
    kmers = {k: 0 for k in all_string_permutations(ALPHABET, KMER_LENGTH)}
    paths = []
    for home_dir, variant_name in zip(VARIANT_HOME_DIRS, VARIANT_NAMES):
        assert os.path.isdir(home_dir), "{} is not a directory".format(home_dir)
        home_dir_paths = os.listdir(home_dir)
        tmp_paths = [home_dir + x + variant_name for x in home_dir_paths if os.path.exists(home_dir + x + variant_name)]
        assert len(tmp_paths) > 0, "Check inputs, there are no paths which exist: {}".format(home_dir)
        paths.append(tmp_paths)

    def get_kmer(chromosome, pos, strand):
        try:
            seq = rh.get_sequence(chromosome, (pos - KMER_LENGTH) + 1, pos + KMER_LENGTH)
            if strand == "-":
                seq = rc.reverse_complement(seq)
            if positions_data is not False:
                replace = read_pos_data[(read_pos_data["contig"] == chromosome) &
                                        (read_pos_data["reference_index"] == pos) &
                                        (read_pos_data["strand"] == strand)]
                if not replace.empty:
                    seq = seq[:KMER_LENGTH - 1] + replace.iloc[0]["replace"] + seq[KMER_LENGTH:]

        except Exception as e:
            print(e, chromosome, pos, strand)
        return seq

    void = '-'
    fill = '#'
    n_spaces = 100
    n_files = 0
    for variant_set, name in zip(paths, NAMES):
        n_paths = len(variant_set)
        count = n_spaces / n_paths
        increaseCount = 0
        print("Starting on {}".format(name))
        local_kmers = {k: 0 for k in all_string_permutations(ALPHABET, KMER_LENGTH)}
        for variant_path in variant_set:
            print('[' + (fill * int(increaseCount)) + (void * int(n_spaces - increaseCount)) + '] ' + str(
                int(increaseCount)) + '%', end='\r')
            increaseCount += count
            variant_data = pd.read_csv(variant_path)
            if positions_data is not False:
                i1 = variant_data.set_index(keys).index
                read_pos_data = positions_data[i2.isin(i1)]
            #         read_id            028a34d4-2a7a-44e7-ab23-305915996ec8
            #         contig                                          RDN18-1
            #         reference_index                                     973
            #         strand                                                +
            #         variants                                             Aa
            #         prob1                                          0.986967
            #         prob2                                          0.013033
            #         prob3                                               NaN
            variant_data['next_base'] = np.vectorize(get_kmer)(variant_data['contig'], variant_data['reference_index'],
                                                               variant_data['strand'])
            large_kmers = set(variant_data['next_base'])
            for l_kmer in large_kmers:
                for i in range(KMER_LENGTH):
                    k = l_kmer[i:KMER_LENGTH + i]
                    if len(k) == KMER_LENGTH:
                        kmers[k] += 1
                        local_kmers[k] += 1

        print('[' + (fill * int(increaseCount)) + (void * int(n_spaces - increaseCount)) + '] ' + str(
            int(increaseCount)) + '%', end='\n')

        total_zeros = 0
        for x, y in local_kmers.items():
            if y == 0:
                total_zeros += 1
        n_files += n_paths
        print("{} Kmers Covered: {}/{}".format(name, len(local_kmers) - total_zeros, len(local_kmers)))
        print("{} Average coverage: {}".format(name, np.average(list(local_kmers.values())) / n_paths))
        with open(os.path.join(OUTPUT_DIR, name + ".tsv"), 'w') as fh:
            print("\n".join(["\t".join([x, str(y / n_paths)]) for x, y in local_kmers.items()]), file=fh)

    total_zeros = 0
    for x, y in kmers.items():
        if y == 0:
            total_zeros += 1
    print("TOTAL Kmers Covered: {}/{}".format(len(kmers) - total_zeros, len(kmers)))
    print("TOTAL Average coverage: {}".format(np.average(list(kmers.values())) / (n_files / 2)))
    with open(os.path.join(OUTPUT_DIR, "total_" + "_".join(NAMES) + ".tsv"), 'w') as fh:
        print("\n".join(["\t".join([x, str(y / n_files)]) for x, y in kmers.items()]), file=fh)


if __name__ == '__main__':
    print(time_it(main)[1])
