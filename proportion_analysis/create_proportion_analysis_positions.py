#!/usr/bin/env python
"""Create a positions file from reference for calling variants. Give assurances on gap size"""
########################################################################
# File: create_proportion_analysis_positions.py
#  executable: create_proportion_analysis_positions.py
#
# Author: Andrew Bailey
# History: 04/19/21 Created
########################################################################

import os
import pandas as pd
import numpy as np
from functools import reduce
from py3helpers.utils import merge_lists
from py3helpers.utils import list_dir, all_string_permutations, time_it
from py3helpers.seq_tools import ReferenceHandler, ReverseComplement
from scipy.stats import norm, invgauss, entropy
# needs pyranges "conda install -c bioconda pyranges"
import pyranges as pr

OUTPUT_DIR = "/home/ubuntu/bisulfite_methylation_analysis/positions"
REFERENCE = "/home/ubuntu/bisulfite_methylation_analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa"
BED_PATHS = "/home/ubuntu/bisulfite_methylation_analysis/bisulfite_data"

REP1_CPG = os.path.join(BED_PATHS, "chr1_ENCFF279HCL.bed")
REP1_CHG = os.path.join(BED_PATHS, "chr1_ENCFF721BJM.bed")
REP1_CHH = os.path.join(BED_PATHS, "chr1_ENCFF448RTC.bed")
REP2_CPG = os.path.join(BED_PATHS, "chr1_ENCFF835NTC.bed")
REP2_CHG = os.path.join(BED_PATHS, "chr1_ENCFF349NNL.bed")
REP2_CHH = os.path.join(BED_PATHS, "chr1_ENCFF038HXQ.bed")

min_coverage = 10
delta = 6
nb_cpu = 10
percents = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
percents = [10, 20, 30, 40, 50, 60, 70, 80, 90]

all_percents = merge_lists([[x-1, x, x+1] for x in percents])

assert os.path.isdir(OUTPUT_DIR), "{} is not a directory".format(OUTPUT_DIR)
assert os.path.exists(REFERENCE), "{} does not exist".format(REFERENCE)
assert os.path.exists(REP1_CPG), "{} does not exist".format(REP1_CPG)
assert os.path.exists(REP1_CHG), "{} does not exist".format(REP1_CHG)
assert os.path.exists(REP1_CHH), "{} does not exist".format(REP1_CHH)
assert os.path.exists(REP2_CPG), "{} does not exist".format(REP2_CPG)
assert os.path.exists(REP2_CHG), "{} does not exist".format(REP2_CHG)
assert os.path.exists(REP2_CHH), "{} does not exist".format(REP2_CHH)

# reference handler and reverse complement handler
rh = ReferenceHandler(REFERENCE)
rc = ReverseComplement()
chromosome_strings = rh.fasta.references[:25]
chromosome_data = {chromosome: rh.get_sequence(chromosome, 0, rh.get_chr_sequence_length(chromosome))
                   for chromosome in rh.fasta.references[:25]}


def get_base(chromosome, pos, strand):
    try:
        base = chromosome_data[chromosome][pos]
        if strand == "-":
            return rc.complement(base)
    except:
        print(chromosome, pos, strand)
    return base


# def get_kmer(pos, strand, offset):
#     start = 0 + offset
#     end = 6 - offset
#     kmer = chr1[pos-start:pos+end]
#     if strand == "-":
#         return rc.complement(kmer)
#     return kmer


##########################################################################################
# Read in Data
##########################################################################################
cpg_rep1 = pd.read_csv(REP1_CPG,
                       names=["Chromosome", "start", "stop", "name", "score", "Strand", "start2", "end2", "color",
                              "coverage", "percentage"],
                       sep="\t",
                       usecols=["Chromosome", "start", "Strand", "coverage", "percentage"])
cpg_rep2 = pd.read_csv(REP2_CPG,
                       names=["Chromosome", "start", "stop", "name", "score", "Strand", "start2", "end2", "color",
                              "coverage", "percentage"],
                       sep="\t",
                       usecols=["Chromosome", "start", "Strand", "coverage", "percentage"])
cpg_data = pd.merge(cpg_rep1, cpg_rep2, on=["Chromosome", "start", "Strand"], how='outer', copy=False,
                    suffixes=("_rep1", "_rep2"))

cpg_rep1 = []
cpg_rep2 = []


chg_rep1 = pd.read_csv(REP1_CHG,
                       names=["Chromosome", "start", "stop", "name", "score", "Strand", "start2", "end2", "color",
                              "coverage", "percentage"],
                       sep="\t",
                       usecols=["Chromosome", "start", "Strand", "coverage", "percentage"])
chg_rep2 = pd.read_csv(REP2_CHG,
                       names=["Chromosome", "start", "stop", "name", "score", "Strand", "start2", "end2", "color",
                              "coverage", "percentage"],
                       sep="\t",
                       usecols=["Chromosome", "start", "Strand", "coverage", "percentage"])
chg_data = pd.merge(chg_rep1, chg_rep2, on=["Chromosome", "start", "Strand"], how='outer', copy=False,
                    suffixes=("_rep1", "_rep2"))
chg_rep1 = []
chg_rep2 = []


chh_rep1 = pd.read_csv(REP1_CHH,
                       names=["Chromosome", "start", "stop", "name", "score", "Strand", "start2", "end2", "color",
                              "coverage", "percentage"],
                       sep="\t",
                       usecols=["Chromosome", "start", "Strand", "coverage", "percentage"])

chh_rep2 = pd.read_csv(REP2_CHH,
                       names=["Chromosome", "start", "stop", "name", "score", "Strand", "start2", "end2", "color",
                              "coverage", "percentage"],
                       sep="\t",
                       usecols=["Chromosome", "start", "Strand", "coverage", "percentage"])

chh_data = pd.merge(chh_rep1, chh_rep2, on=["Chromosome", "start", "Strand"], how='outer', copy=False,
                    suffixes=("_rep1", "_rep2"))
chh_rep1 = []
chh_rep2 = []

all_data = pd.concat([cpg_data, chg_data, chh_data], ignore_index=True)
cpg_data = []
chg_data = []
chh_data = []

possible_mods = all_data[((all_data.coverage_rep1 >= min_coverage) &
                          (all_data.coverage_rep2 >= min_coverage) &
                          (all_data.percentage_rep1.isin(all_percents)) &
                          (all_data.percentage_rep1 == all_data.percentage_rep2))]

possible_mods['Start'] = possible_mods['start'] - delta
possible_mods['End'] = possible_mods['start'] + delta
possible_mods['find'] = "C"
possible_mods['replace'] = "M"
possible_mods['ambig'] = "Y"

# output_path = os.path.join(BED_PATHS, "filtered_chr1_all.bed")
# possible_mods.to_csv(output_path, sep="\t", index=False)

for x in percents:
    output_path = os.path.join(BED_PATHS, f"filtered_chr1_{x}.bed")
    temp_mods = possible_mods[possible_mods["percentage_rep2"].isin([x-1, x, x+1])]
    temp_mods.to_csv(output_path, sep="\t", index=False)
    output_path = os.path.join(OUTPUT_DIR, f"filtered_chr1_{x}.positions")
    temp_mods.to_csv(output_path, sep="\t", index=False,
                             columns=["Chromosome", "start", "Strand", "find", "replace"], header=False)
    output_path = os.path.join(OUTPUT_DIR, f"filtered_chr1_{x}_variant.positions")
    temp_mods.to_csv(output_path, sep="\t", index=False,
                     columns=["Chromosome", "start", "Strand", "find", "ambig"], header=False)

##########################################################################################
