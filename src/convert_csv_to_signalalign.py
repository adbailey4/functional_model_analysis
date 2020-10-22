#!/usr/bin/env python
"""Convert csv model to signalalign model"""
########################################################################
# File: convert_csv_to_signalalign.py
#  executable: convert_csv_to_signalalign.py
#
# Author: Andrew Bailey
# History: 05/04/20 Created
########################################################################

import os
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from signalalign.hiddenMarkovModel import HmmModel
from py3helpers.utils import list_dir, all_string_permutations, time_it
from py3helpers.multiprocess import BasicService, run_service


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--dir', '-d', action='store',
                        dest='dir', required=True, type=str, default=None,
                        help="Path to directory of csv files to convert to signalalign models")
    parser.add_argument('--output_dir', '-o', action='store',
                        dest='output_dir', required=True, type=str, default=None,
                        help="Path to directory to output signalalign models")
    parser.add_argument('--base_model', '-m', action='store',
                        dest='base_model', required=True, type=str, default=None,
                        help="Path to base signalalign model")
    parser.add_argument('--num_threads', '-t', action='store',
                        dest='num_threads', required=False, type=int, default=4,
                        help="Number of threads to run")
    parser.add_argument('--rna', '-r', action='store_true',
                        dest='rna', required=False,
                        help="Boolean option if models are for rna")

    args = parser.parse_args()
    return args


def convert_csv_to_sa_model(csv_file, output_dir, transition_probs, state_number=3):
    output_path = os.path.join(output_dir, os.path.splitext(os.path.basename(csv_file))[0]+".model")
    data = pd.read_csv(csv_file, names=["kmer", "mean", "sd"])
    alphabet = "".join(sorted(set("".join(data["kmer"]))))
    kmer_length = len(data["kmer"].iloc[0])
    alphabet_size = len(alphabet)
    new_kmers = all_string_permutations(alphabet, length=kmer_length)
    with open(output_path, 'w') as f:

        # line 0
        f.write("{stateNumber}\t{alphabetSize}\t{alphabet}\t{kmerLength}\n"
                "".format(stateNumber=state_number, alphabetSize=alphabet_size,
                          alphabet=alphabet, kmerLength=kmer_length))
        # line 1 transitions
        for i in range(state_number * state_number):
            f.write("{transition}\t".format(transition=str(transition_probs[i])))
        # likelihood
        f.write("{}\n".format(str(0)))

        # line 2 Event Model
        for kmer in new_kmers:
            # k_index = HmmModel._get_kmer_index(kmer, alphabet, kmer_length, alphabet_size)
            kmer_data = data[data["kmer"] == kmer]
            assert kmer == kmer_data["kmer"].iloc[0], "The input csv model is not sorted or something else is very wrong. Check inputs"
            f.write("{level_mean}\t{level_sd}\t{noise_mean}\t{noise_sd}\t{noise_lambda}\t"
                    "".format(level_mean=kmer_data["mean"].iloc[0],
                              level_sd=kmer_data["sd"].iloc[0],
                              noise_mean=0,
                              noise_sd=0,
                              noise_lambda=0))
        f.write("\n")
    return output_path


def main():
    args = parse_args()
    print(args)
    assert os.path.isdir(args.dir), "{} is not a directory".format(args.dir)
    assert os.path.isdir(args.output_dir), "{} is not a directory".format(args.output_dir)
    assert os.path.exists(args.base_model), "{} does not exist".format(args.base_model)
    csv_files = list_dir(args.dir, ext="csv")
    model = HmmModel(args.base_model, rna=args.rna)
    transition_probs = model.transitions

    extra_args = {"output_dir": args.output_dir,
                  "transition_probs": transition_probs,
                  "state_number": 3}
    service = BasicService(convert_csv_to_sa_model, service_name="multiprocess_convert_csv_to_sa_model")
    total, failure, messages, output = run_service(service.run, csv_files,
                                                   extra_args, ["csv_file"], args.num_threads)


if __name__ == '__main__':
    print(time_it(main)[1])
