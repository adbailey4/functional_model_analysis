#!/usr/bin/env python
"""Get the minimum number of reads that cover all kmers given a positions file"""
########################################################################
# File: minimum_kmer_covering_reads.py
#  executable: minimum_kmer_covering_reads.py
#
# Author: Andrew Bailey
# History: 05/24/20 Created
########################################################################

import os
import sys
import pandas as pd
import numpy as np
from argparse import ArgumentParser
import pysam
from contextlib import closing

from py3helpers.utils import all_string_permutations, merge_lists, time_it
from py3helpers.seq_tools import ReferenceHandler, ReverseComplement
from py3helpers.multiprocess import run_service, BasicService


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--bam', '-b', action='store',
                        dest='bam', required=True, type=str, default=None,
                        help="Path to bam file")
    parser.add_argument('--positions_file', '-p', action='store',
                        dest='positions_file', required=True, type=str, default=None,
                        help="Path to positions_file")
    parser.add_argument('--output_dir', '-o', action='store',
                        dest='output_dir', required=True, type=str, default=None,
                        help="Path to output directory")
    parser.add_argument('--alphabet', '-a', action='store',
                        dest='alphabet', required=False, type=str, default="ATGC",
                        help="Set of characters to cover")
    parser.add_argument('--kmer_length', '-l', action='store',
                        dest='kmer_length', required=False, type=int, default=6,
                        help="Length of the kmer to cover")
    parser.add_argument('--reference', '-r', action='store',
                        dest='reference', required=False, type=str, default=None,
                        help="Reference if MD flag is not set in BAM")

    args = parser.parse_args()
    return args


class Kmer(object):
    def __init__(self, kmer):
        self.kmer = kmer
        self.read_indices = []

    def add_read_index(self, read_index):
        self.read_indices.append(read_index)

    def remove_read_index(self, read_index):
        for i, x in enumerate(self.read_indices):
            if x == read_index:
                self.read_indices.pop(i)
                break

    def get_read_index(self):
        if len(self.read_indices) == 0:
            return None
        else:
            return self.read_indices[0]


class Read(object):
    def __init__(self, read_id):
        self.read_id = read_id
        self.kmers = []

    def add_kmer(self, kmer):
        self.kmers.append(kmer)


class KmerMap(object):
    def __init__(self, alphabet="ACGT", length=6):
        self.alphabet = "".join(sorted(alphabet))
        self.length = length
        self.kmers = [x for x in all_string_permutations(self.alphabet, self.length)]
        self.reads = []
        self.kmer_classes = {x: Kmer(x) for x in self.kmers}
        self.kmer_counts = {x: 0 for x in self.kmers}

    def add_read(self, read):
        read_index = len(self.reads)
        self.reads.append(read)
        for kmer in read.kmers:
            self.kmer_classes[kmer].add_read_index(read_index)
            self.kmer_counts[kmer] += 1

    def get_min_kmer(self):
        return min(self.kmer_counts, key=self.kmer_counts.get)

    def get_non_zero_min_kmer(self):
        min_count = np.inf
        return_kmer = None
        for kmer, count in self.kmer_counts.items():
            if min_count > count > 0:
                min_count = count
                return_kmer = kmer
        return return_kmer

    def get_non_zero_min_kmer_in_kmers(self, kmers):
        min_count = np.inf
        return_kmer = None
        for kmer, count in self.kmer_counts.items():
            if kmer in kmers:
                if min_count > count > 0:
                    min_count = count
                    return_kmer = kmer
        return return_kmer

    def get_zero_kmers(self):
        zero_kmers = []
        for kmer, count in self.kmer_counts.items():
            if count == 0:
                zero_kmers.append(kmer)
        return zero_kmers

    def remove_read(self, read_index):
        read = self.reads[read_index]
        for kmer in read.kmers:
            self.kmer_classes[kmer].remove_read_index(read_index)
            self.kmer_counts[kmer] -= 1

    def get_threshold_covered_kmers(self, threshold=1):
        kmers = []
        for kmer, count in self.kmer_counts.items():
            if count >= threshold:
                kmers.append(kmer)
        return kmers

    def get_threshold_uncovered_kmers(self, threshold=1):
        kmers = []
        for kmer, count in self.kmer_counts.items():
            if count < threshold:
                kmers.append(kmer)
        return kmers

    def get_read(self, kmer):
        read_index = self.kmer_classes[kmer].get_read_index()
        if read_index is None:
            return read_index, None
        else:
            return read_index, self.reads[read_index]


def main():
    args = parse_args()
    assert os.path.isdir(args.output_dir), "{} is not a directory".format(args.output_dir)
    assert os.path.exists(args.bam), "{} does not exist".format(args.bam)
    assert os.path.exists(args.positions_file), "{} does not exist".format(args.positions_file)

    fasta_handle = None
    if args.reference is not None:
        assert os.path.exists(args.reference), "{} does not exist".format(args.reference)
        fasta_handle = ReferenceHandler(args.reference)

    rc = ReverseComplement()
    positions_data = pd.read_csv(args.positions_file, names=["chr", "start", "strand", "find", "replace"], sep="\t")
    km = KmerMap(args.alphabet, args.kmer_length)
    counter = 0

    def get_kmer(sequence, pos, start_pos, strand):
        try:
            base = sequence[(pos - (args.kmer_length-1)) - start_pos:(pos + args.kmer_length) - start_pos]
            if strand == "-":
                return rc.complement(base)
            return base
        except Exception as e:
            print(e, sequence, pos, start_pos)

    # def get_ref_base(chromosome, start_pos, strand):
    #     try:
    #         base = fasta_handle.get_sequence(chromosome_name=chromosome, start=start_pos, stop=start_pos + 1)
    #         if strand == "-":
    #             return rc.complement(base)
    #         return base
    #     except Exception as e:
    #         print(e, fasta_handle, chromosome, start_pos, strand)
    #
    # def get_base(sequence, pos, start_pos, reversed):
    #     try:
    #         base = sequence[pos - start_pos]
    #         if reversed:
    #             return rc.complement(base)
    #         return base
    #     except Exception as e:
    #         print(e, sequence, pos, start_pos)

    def get_covered_kmers(positions_data1, read_name1, ref_sequence1, ref_name1, strand1, ref_start1, ref_end1):
        this_positions_data = positions_data1.loc[(positions_data1["chr"] == ref_name1) &
                                                  (positions_data1["strand"] == strand1) &
                                                  (positions_data1["start"] >= ref_start1) &
                                                  (positions_data1["start"] <= ref_end1)]
        if this_positions_data.empty:
            return None
        kmer_lists = np.vectorize(get_kmer)(ref_sequence1,
                                            this_positions_data['start'],
                                            ref_start1,
                                            strand1)
        kmer_subset_lists1 = merge_lists([[kmer[i:i + args.kmer_length] for i in range(args.kmer_length) if len(kmer[i:i + args.kmer_length]) == args.kmer_length and set(kmer[i:i + args.kmer_length]) <= set(args.alphabet)] for kmer in kmer_lists])
        return read_name1, kmer_subset_lists1

    def meta_get_covered_kmers(positions, all_args1):
        data_to_return = []
        for args1 in all_args1:
            data = get_covered_kmers(positions, *args1)
            if data is not None:
                data_to_return.append(data)
        return data_to_return

    all_args = []
    with closing(pysam.AlignmentFile(args.bam, 'rb' if args.bam.endswith("bam") else 'r')) as aln:
        for aligned_segment in aln.fetch(until_eof=True):
            if counter > 2000:
                break
            try:
                if not aligned_segment.has_tag('MD'):
                    if fasta_handle is None:
                        raise Exception("Need to specify --reference if MD flag is not set")
                    else:
                        ref_sequence = fasta_handle.get_sequence(chromosome_name=aligned_segment.reference_name,
                                                                 start=aligned_segment.reference_start,
                                                                 stop=aligned_segment.reference_end)
                else:
                    ref_sequence = aligned_segment.get_reference_sequence().upper()

                read_name = aligned_segment.qname.split("_")[0]
                ref_name = aligned_segment.reference_name
                ref_start = aligned_segment.reference_start
                ref_end = aligned_segment.reference_end
                reversed_read = aligned_segment.is_reverse
                if reversed_read:
                    strand = "-"
                else:
                    strand = "+"

                all_args.append([read_name, ref_sequence, ref_name, strand, ref_start, ref_end])

                counter += 1

            except Exception as e:
                print(e, file=sys.stderr)

    n_processes = 4
    print("starting on {} reads".format(len(all_args)))
    list_of_args = [all_args[x::n_processes] for x in range(n_processes)]
    extra_args = {"positions": positions_data}
    service = BasicService(meta_get_covered_kmers, service_name="multiprocess_meta_get_covered_kmers")
    total, failure, messages, output = run_service(service.run,
                                                   list_of_args,
                                                   extra_args,
                                                   ["all_args1"],
                                                   n_processes)
    # print(pd.concat(output, ignore_index=True))

    all_data = merge_lists(output)
    print("number of reads: ", len(all_data))
    for read_name, kmer_subset_lists in all_data:
        # print(read_name, kmer_subset_lists)
        r = Read(read_name)
        for kmer in kmer_subset_lists:
            r.add_kmer(kmer)
        km.add_read(r)

    keep_kmer_map = KmerMap(args.alphabet, args.kmer_length)

    print("number of zero covered kmers: ", len(km.get_zero_kmers()))
    curr_threshold = 1
    find_kmers = keep_kmer_map.get_threshold_uncovered_kmers(threshold=curr_threshold)
    iteration = 1
    increase_threshold = True
    while increase_threshold:
        while len(find_kmers) > 0:
            print(iteration, len(find_kmers))
            next_kmer = km.get_non_zero_min_kmer_in_kmers(find_kmers)
            if next_kmer is None:
                print("No more reads to cover found kmers: threshold {}".format(curr_threshold))
                increase_threshold = False
                break
            next_read_index, next_read = km.get_read(next_kmer)
            if next_read is None:
                print("Whoops, something is wrong")
                break
            keep_kmer_map.add_read(next_read)
            km.remove_read(next_read_index)
            find_kmers = keep_kmer_map.get_threshold_uncovered_kmers(threshold=curr_threshold)
            iteration += 1
        if len(find_kmers) == 0:
            print("Found reads covering all kmers at threshold {}".format(curr_threshold))
        file_path = os.path.join(args.output_dir, "{}_reads_covering_kmers_with_threshold_{}.txt".format("all" if increase_threshold else "some", curr_threshold))
        with open(file_path, "w") as fh:
            print("\n".join([read.read_id for read in keep_kmer_map.reads]), file=fh)
        kmer_counts_file_path = os.path.join(args.output_dir, "{}_kmer_counts_with_threshold_{}.txt".format("all" if increase_threshold else "some", curr_threshold))
        with open(kmer_counts_file_path, "w") as fh:
            print("\n".join(["\t".join([kmer, str(count)]) for kmer, count in keep_kmer_map.kmer_counts.items()]), file=fh)

    curr_threshold += 1


if __name__ == '__main__':
    print(time_it(main)[1])
