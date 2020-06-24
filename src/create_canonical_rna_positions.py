#!/usr/bin/env python
"""Create a positions file from reference for calling variants. Give assurances on gap size"""
########################################################################
# File: create_canonical_positions.py
#  executable: create_canonical_positions.py
#
# Author: Andrew Bailey
# History: 05/11/20 Created
########################################################################

import os
import pandas as pd
import numpy as np
from functools import reduce

from py3helpers.utils import list_dir, all_string_permutations, time_it
from py3helpers.seq_tools import ReferenceHandler, ReverseComplement
from scipy.stats import norm, invgauss, entropy
# needs pyranges "conda install -c bioconda pyranges"
import pyranges as pr

OUTPUT_DIR = "/home/ubuntu/bisulfite_methylation_analysis/positions"
REFERENCE = "/home/ubuntu/bisulfite_methylation_analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
p_lambda = 50
delta = 6

assert os.path.isdir(OUTPUT_DIR), "{} is not a directory".format(OUTPUT_DIR)
assert os.path.exists(REFERENCE), "{} does not exist".format(REFERENCE)

# reference handler and reverse complement handler
rh = ReferenceHandler(REFERENCE)
rc = ReverseComplement()
transcript_strings = rh.fasta.references
transcript_data = {transcript: rh.get_sequence(transcript, 0, rh.get_chr_sequence_length(transcript))
                   for transcript in rh.fasta.references}


def get_base(transcript, pos):
    try:
        base = transcript_data[transcript][pos]
    except:
        print(chromosome, pos)
    return base


# def get_kmer(pos, strand, offset):
#     start = 0 + offset
#     end = 6 - offset
#     kmer = chr1[pos-start:pos+end]
#     if strand == "-":
#         return rc.complement(kmer)
#     return kmer


##########################################################################################
# confident canonical calls
##########################################################################################


transcripts = []
strands = []
starts = []
ends = []

# randomly create positions to call on chromosomes
for transcript in transcript_strings:
    chr_len = rh.get_chr_sequence_length(transcript)
    s = np.random.poisson(p_lambda, int(chr_len / p_lambda))
    positions = np.cumsum(s)
    indexes_larger_than_ref, = np.nonzero(positions >= chr_len)
    # print(type(indexes_larger_than_ref))
    if len(indexes_larger_than_ref) > 0:
        max_index = np.min(indexes_larger_than_ref)
        positions = positions[:max_index]
    start = positions - delta
    end = positions + delta
    transcripts.extend([transcript] * len(positions))
    starts.extend(list(start))
    ends.extend(list(end))
    print(transcript)

# random positions as pyranges
canonical_positions = pr.from_dict({"Chromosome": transcripts, "Start": starts, "End": ends})

# remove cut off canonical ranges
canonical_positions = canonical_positions[(canonical_positions.End - canonical_positions.Start) == (delta * 2)]
canonical_positions = canonical_positions.df
# recreate original position for left over data
canonical_positions["midpoint"] = canonical_positions.Start + delta
# get bases for each position (all should be C but this is worth a double check)
canonical_positions['find'] = np.vectorize(get_base)(canonical_positions['Chromosome'], canonical_positions['midpoint'], canonical_positions['Strand'])
canonical_positions['replace'] = canonical_positions['find']
canonical_positions['ambig'] = "X"
canonical_positions = canonical_positions[canonical_positions["find"] != "N"]
# output data
output_path = os.path.join(OUTPUT_DIR, "transcript_canonical.positions")
ambig_output_path = os.path.join(OUTPUT_DIR, "transcript_canonical_variant.positions")
canonical_positions.to_csv(output_path, sep="\t", index=False,
                           columns=["Chromosome", "midpoint", "Strand", "find", "replace"], header=False)
canonical_positions.to_csv(ambig_output_path, sep="\t", index=False,
                           columns=["Chromosome", "midpoint", "Strand", "find", "ambig"], header=False)

