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
REP1_CPG = "/home/ubuntu/bisulfite_methylation_analysis/bisulfite_data/ENCFF279HCL.bed"
REP1_CHG = "/home/ubuntu/bisulfite_methylation_analysis/bisulfite_data/ENCFF721BJM.bed"
REP1_CHH = "/home/ubuntu/bisulfite_methylation_analysis/bisulfite_data/ENCFF448RTC.bed"
REP2_CPG = "/home/ubuntu/bisulfite_methylation_analysis/bisulfite_data/ENCFF835NTC.bed"
REP2_CHG = "/home/ubuntu/bisulfite_methylation_analysis/bisulfite_data/ENCFF349NNL.bed"
REP2_CHH = "/home/ubuntu/bisulfite_methylation_analysis/bisulfite_data/ENCFF038HXQ.bed"

min_coverage = 10
delta = 6
nb_cpu = 10

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
chg_rep1 = pd.read_csv(REP1_CHG,
                       names=["Chromosome", "start", "stop", "name", "score", "Strand", "start2", "end2", "color",
                              "coverage", "percentage"],
                       sep="\t",
                       usecols=["Chromosome", "start", "Strand", "coverage", "percentage"])
chh_rep1 = pd.read_csv(REP1_CHH,
                       names=["Chromosome", "start", "stop", "name", "score", "Strand", "start2", "end2", "color",
                              "coverage", "percentage"],
                       sep="\t",
                       usecols=["Chromosome", "start", "Strand", "coverage", "percentage"])
cpg_rep2 = pd.read_csv(REP2_CPG,
                       names=["Chromosome", "start", "stop", "name", "score", "Strand", "start2", "end2", "color",
                              "coverage", "percentage"],
                       sep="\t",
                       usecols=["Chromosome", "start", "Strand", "coverage", "percentage"])
chg_rep2 = pd.read_csv(REP2_CHG,
                       names=["Chromosome", "start", "stop", "name", "score", "Strand", "start2", "end2", "color",
                              "coverage", "percentage"],
                       sep="\t",
                       usecols=["Chromosome", "start", "Strand", "coverage", "percentage"])
chh_rep2 = pd.read_csv(REP2_CHH,
                       names=["Chromosome", "start", "stop", "name", "score", "Strand", "start2", "end2", "color",
                              "coverage", "percentage"],
                       sep="\t",
                       usecols=["Chromosome", "start", "Strand", "coverage", "percentage"])

# process cpg confident positions
cpg_data = pd.merge(cpg_rep1, cpg_rep2, on=["Chromosome", "start", "Strand"], how='outer', copy=False,
                    suffixes=("_rep1", "_rep2"))
chg_data = pd.merge(chg_rep1, chg_rep2, on=["Chromosome", "start", "Strand"], how='outer', copy=False,
                    suffixes=("_rep1", "_rep2"))
chh_data = pd.merge(chh_rep1, chh_rep2, on=["Chromosome", "start", "Strand"], how='outer', copy=False,
                    suffixes=("_rep1", "_rep2"))
cpg_data["cpg"] = True
chg_data["cpg"] = False
chh_data["cpg"] = False

cpg_rep1 = []
chg_rep1 = []
chh_rep1 = []
cpg_rep2 = []
chg_rep2 = []
chh_rep2 = []

all_data = pd.concat([cpg_data, chg_data, chh_data], ignore_index=True)
cpg_data = []
chg_data = []
chh_data = []
all_data["confident_unmod"] = ((all_data.coverage_rep1 >= min_coverage) & (all_data.coverage_rep2 >= min_coverage) & (
        all_data.percentage_rep1 <= max_canonical_percentage) & (
                                       all_data.percentage_rep2 <= max_canonical_percentage))
all_data["confident_mod"] = ((all_data.coverage_rep1 >= min_coverage) & (all_data.coverage_rep2 >= min_coverage) & (
        all_data.percentage_rep1 >= min_mod_percentage) & (all_data.percentage_rep2 >= min_mod_percentage))
all_data["confident_cpg_mod"] = all_data["confident_mod"] & all_data["cpg"]
all_data['Start'] = all_data['start'] - delta
all_data['End'] = all_data['start'] + delta
all_data['find'] = "C"

##########################################################################################
# confident mC mods
##########################################################################################
# get interval data for non mod positions
excluded_intervals = pr.PyRanges(all_data[~all_data.confident_unmod & ~all_data.confident_mod]).merge(strand=True)
included_intervals = pr.PyRanges(all_data[all_data["confident_mod"]])
all_mod_intervals = included_intervals.subtract(excluded_intervals, nb_cpu=nb_cpu).df
# isolated mods
all_mod_intervals = all_mod_intervals[(all_mod_intervals.End - all_mod_intervals.Start) == (delta * 2)]
all_mod_intervals["midpoint"] = all_mod_intervals.Start + delta
# get bases for each position (all should be C but this is worth a double check)
all_mod_intervals['find'] = np.vectorize(get_base)(all_mod_intervals['Chromosome'], all_mod_intervals['midpoint'],
                                                   all_mod_intervals['Strand'])
all_mod_intervals['replace'] = "M"
all_mod_intervals['ambig'] = "Y"
all_mod_intervals = all_mod_intervals[all_mod_intervals["find"] != "N"]
output_path = os.path.join(OUTPUT_DIR, "all_mC.positions")
ambig_output_path = os.path.join(OUTPUT_DIR, "all_mC_variant.positions")
all_mod_intervals.to_csv(output_path, sep="\t", index=False,
                         columns=["Chromosome", "midpoint", "Strand", "find", "replace"], header=False)
all_mod_intervals.to_csv(ambig_output_path, sep="\t", index=False,
                         columns=["Chromosome", "midpoint", "Strand", "find", "ambig"], header=False)
# just cpg
output_path = os.path.join(OUTPUT_DIR, "mc_cpg.positions")
ambig_output_path = os.path.join(OUTPUT_DIR, "mc_cpg_variant.positions")
all_mod_intervals[all_mod_intervals["cpg"]].to_csv(output_path, sep="\t", index=False,
                                                   columns=["Chromosome", "midpoint", "Strand", "find", "replace"],
                                                   header=False)
all_mod_intervals[all_mod_intervals["cpg"]].to_csv(ambig_output_path, sep="\t", index=False,
                                                   columns=["Chromosome", "midpoint", "Strand", "find", "ambig"],
                                                   header=False)
# just not cpg
output_path = os.path.join(OUTPUT_DIR, "mc_not_cpg.positions")
ambig_output_path = os.path.join(OUTPUT_DIR, "mc_not_cpg_variant.positions")
all_mod_intervals[~all_mod_intervals["cpg"]].to_csv(output_path, sep="\t", index=False,
                                                    columns=["Chromosome", "midpoint", "Strand", "find", "replace"],
                                                    header=False)
all_mod_intervals[~all_mod_intervals["cpg"]].to_csv(ambig_output_path, sep="\t", index=False,
                                                    columns=["Chromosome", "midpoint", "Strand", "find", "ambig"],
                                                    header=False)
# clean up mods
all_mod_intervals = []
excluded_intervals = []
included_intervals = []

##########################################################################################
# confident cpg canonical mods
##########################################################################################
# create range from non confident unmodified positions
possible_mod_intervals = pr.PyRanges(all_data[~all_data["confident_unmod"]]).merge(strand=True)
included_intervals = pr.PyRanges(all_data[all_data["confident_unmod"]])
all_cpg_unmod_intervals_pr = included_intervals.subtract(possible_mod_intervals, nb_cpu=nb_cpu)
all_cpg_unmod_intervals = all_cpg_unmod_intervals_pr.df  # 427,387,017
all_cpg_unmod_intervals = all_cpg_unmod_intervals[(all_cpg_unmod_intervals.End - all_cpg_unmod_intervals.Start) == (delta * 2)]
all_cpg_unmod_intervals_pr = pr.PyRanges(all_cpg_unmod_intervals)
all_cpg_unmod_intervals_pr = all_cpg_unmod_intervals_pr.merge(strand=True)
all_cpg_unmod_intervals = all_cpg_unmod_intervals_pr.df
all_cpg_unmod_intervals = all_cpg_unmod_intervals[(all_cpg_unmod_intervals.End - all_cpg_unmod_intervals.Start) == (delta * 2)]

# recreate original position for left over data
all_cpg_unmod_intervals["midpoint"] = all_cpg_unmod_intervals.Start + delta
# get bases for each position (all should be C but this is worth a double check)
all_cpg_unmod_intervals['find'] = np.vectorize(get_base)(all_cpg_unmod_intervals['Chromosome'], all_cpg_unmod_intervals['midpoint'],
                                                         all_cpg_unmod_intervals['Strand'])
all_cpg_unmod_intervals['replace'] = all_cpg_unmod_intervals['find']
all_cpg_unmod_intervals['ambig'] = "Y"
all_cpg_unmod_intervals = all_cpg_unmod_intervals[all_cpg_unmod_intervals["find"] != "N"]
# output data
output_path = os.path.join(OUTPUT_DIR, "cxx_canonical.positions")
ambig_output_path = os.path.join(OUTPUT_DIR, "cxx_canonical_variant.positions")
all_cpg_unmod_intervals.to_csv(output_path, sep="\t", index=False,
                               columns=["Chromosome", "midpoint", "Strand", "find", "replace"], header=False)
all_cpg_unmod_intervals.to_csv(ambig_output_path, sep="\t", index=False,
                               columns=["Chromosome", "midpoint", "Strand", "find", "ambig"], header=False)

##########################################################################################
# confident canonical mods
##########################################################################################

chromosomes = []
strands = []
starts = []
ends = []

# randomly create positions to call on chromosomes
for chromosome in chromosome_strings:
    chr_len = rh.get_chr_sequence_length(chromosome)
    s = np.random.poisson(p_lambda, int(chr_len / p_lambda))
    positions = np.cumsum(s)
    indexes_larger_than_ref, = np.nonzero(positions >= chr_len)
    # print(type(indexes_larger_than_ref))
    if len(indexes_larger_than_ref) > 0:
        max_index = np.min(indexes_larger_than_ref)
        positions = positions[:max_index]
    start = positions - delta
    end = positions + delta
    chromosomes.extend([chromosome] * (2 * len(positions)))
    strands.extend((["+"] * len(positions)) + (["-"] * len(positions)))
    starts.extend(list(np.append(start, start)))
    ends.extend(list(np.append(end, end)))
    print(chromosome)

# random positions as pyranges
gr2 = pr.from_dict({"Chromosome": chromosomes, "Strand": strands, "Start": starts, "End": ends})

all_cpg_unmod_intervals_pr = pr.PyRanges(all_cpg_unmod_intervals)
# remove unconfident ranges from randomly created positions
canonical_positions = gr2.subtract(possible_mod_intervals)
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
output_path = os.path.join(OUTPUT_DIR, "canonical_minimum_cxx.positions")
ambig_output_path = os.path.join(OUTPUT_DIR, "canonical_minimum_cxx_variant.positions")
canonical_positions.to_csv(output_path, sep="\t", index=False,
                           columns=["Chromosome", "midpoint", "Strand", "find", "replace"], header=False)
canonical_positions.to_csv(ambig_output_path, sep="\t", index=False,
                           columns=["Chromosome", "midpoint", "Strand", "find", "ambig"], header=False)


##########################################################################################
# confident canonical mods plus cpgs
##########################################################################################
canonical_positions = gr2.subtract(possible_mod_intervals)
# remove cut off canonical ranges
canonical_positions = canonical_positions[(canonical_positions.End - canonical_positions.Start) == (delta * 2)]
canonical_positions = canonical_positions.df
canonical_positions = pr.PyRanges(canonical_positions).subtract(pr.PyRanges(all_cpg_unmod_intervals_pr.df.iloc[::6, :]))
canonical_positions = canonical_positions[(canonical_positions.End - canonical_positions.Start) == (delta * 2)]
canonical_positions = canonical_positions.df
canonical_positions["midpoint"] = canonical_positions.Start + delta
# get bases for each position (all should be C but this is worth a double check)
canonical_positions['find'] = np.vectorize(get_base)(canonical_positions['Chromosome'], canonical_positions['midpoint'],
                                                     canonical_positions['Strand'])
canonical_positions = canonical_positions[canonical_positions["find"] != "N"]
canonical_positions['replace'] = canonical_positions['find']
canonical_positions['ambig'] = "X"

# canonical_positions = gr2.subtract(possible_mod_intervals).df
canonical_positions = pd.concat([canonical_positions, all_cpg_unmod_intervals_pr.df.iloc[::6, :]], ignore_index=True)
canonical_positions['ambig'] = "X"

# remove cut off canonical ranges
# recreate original position for left over data
# output data
output_path = os.path.join(OUTPUT_DIR, "canonical_added_cxx.positions")
ambig_output_path = os.path.join(OUTPUT_DIR, "canonical_added_cxx_variant.positions")
canonical_positions.to_csv(output_path, sep="\t", index=False,
                           columns=["Chromosome", "midpoint", "Strand", "find", "replace"], header=False)
canonical_positions.to_csv(ambig_output_path, sep="\t", index=False,
                           columns=["Chromosome", "midpoint", "Strand", "find", "ambig"], header=False)


