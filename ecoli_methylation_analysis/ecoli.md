## Ecoli analysis
PLAN: Download Fast5 files, call with Megalodon, use megalodon calls to label positions

* Data from `https://www.nature.com/articles/s41592-021-01109-3`
* Raw Fast5 data `https://sra-download.ncbi.nlm.nih.gov/traces/sra59/SRZ/010032/SRR10032546/MinION_EC_NAT.tar.gz`

* Download data
```
wget https://sra-download.ncbi.nlm.nih.gov/traces/sra59/SRZ/010032/SRR10032546/MinION_EC_NAT.tar.gz
tar -xzf MinION_EC_NAT.tar.gz
# reference from https://www.ncbi.nlm.nih.gov/nuccore/U00096

nanopolish extract /home/ubuntu/ecoli_methylation_analysis/MinION_EC_NAT --recurse --fastq --output /home/ubuntu/ecoli_methylation_analysis/fastq/MinION_EC_NAT.fastq
minimap2 --MD -t 18 -ax map-ont /home/ubuntu/ecoli_methylation_analysis/reference/ecoli.fa /home/ubuntu/ecoli_methylation_analysis/fastq/MinION_EC_NAT.fastq | samtools view -@ 10 -bS - | samtools sort -@ 10 - > /home/ubuntu/ecoli_methylation_analysis/fastq/MinION_EC_NAT.sorted.bam && samtools view -@ 40 -bSF 2308 /home/ubuntu/ecoli_methylation_analysis/fastq/MinION_EC_NAT.sorted.bam > /home/ubuntu/ecoli_methylation_analysis/fastq/MinION_EC_NAT.2308.sorted.bam && samtools index /home/ubuntu/ecoli_methylation_analysis/fastq/MinION_EC_NAT.2308.sorted.bam
samtools coverage -A -w 50 -r U00096.3 /home/ubuntu/ecoli_methylation_analysis/fastq/MinION_EC_NAT.2308.sorted.bam
# remove_sa_analyses -d /home/ubuntu/ecoli_methylation_analysis/MinION_EC_NAT --analysis --threads 18
```
![](ecoli_coverage.png)

* Run megalodon as our labels
```
# megalodon on p3.2xlarge
megalodon /home/ubuntu/ecoli_methylation_analysis/MinION_EC_NAT \
--outputs mod_mappings mods \
--reference /home/ubuntu/ecoli_methylation_analysis/reference/ecoli.fa \
--processes 8 --overwrite --guppy-server-path guppy_basecall_server \
--output-directory /home/ubuntu/ecoli_methylation_analysis/output \
--guppy-timeout 1000 --guppy-params "--device cuda:0"
```

* Position Filtering

```
import pandas as pd
import os

N = 20
megalodon = "/home/ubuntu/ecoli_methylation_analysis/output/modified_bases.5mC.bed"
BED_PATHS = "/home/ubuntu/ecoli_methylation_analysis/subset_reference/bed_files"
OUTPUT_DIR = "/home/ubuntu/ecoli_methylation_analysis/subset_reference/positions_files"
top_n_dir = "/home/ubuntu/ecoli_methylation_analysis/subset_reference/filtered_positions"

possible_mods = pd.read_csv(megalodon,
                       names=["Chromosome", "start", "stop", "name", "score", "Strand", "start2", "end2", "color",
                              "coverage", "percentage"],
                       sep="\t",
                       usecols=["Chromosome", "start", "Strand", "coverage", "percentage"])
percents = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]


possible_mods['find'] = "C"
possible_mods['replace'] = "M"
possible_mods['ambig'] = "Y"
possible_mods = possible_mods[possible_mods["start"] < 100000]
for x in percents:
    output_path = os.path.join(BED_PATHS, f"filtered_chr1_{x}.bed")
    temp_mods = possible_mods[(possible_mods["percentage"] > x-1) & (possible_mods["percentage"] < x+1)]
    temp_mods.to_csv(output_path, sep="\t", index=False)
    output_path = os.path.join(OUTPUT_DIR, f"filtered_{x}_variant.positions")
    temp_mods.to_csv(output_path, sep="\t", index=False,
                     columns=["Chromosome", "start", "Strand", "find", "ambig"], header=False)
    sample_number = len(temp_mods)
    if N < len(temp_mods):
        sample_number = N 
    sample = temp_mods.sample(sample_number)
    output_path = os.path.join(top_n_dir, f"top_{N}_covered_filtered_{x}.positions")
    sample.to_csv(output_path, sep="\t", index=False,
                     columns=["Chromosome", "start", "Strand", "find", "ambig"], header=False)

```

* Prep data for S3 / kubernetes
  `bash prep_data.sh`

* Run K8 job
  `make ecoli`

* Kmer counts
```
# uncomment ECOLI MOD
python /home/ubuntu/functional_model_analysis/src/coverage_from_variant_calls.py
```


### Kmer analysis
```
import pandas as pd
import os

megalodon = "/home/ubuntu/ecoli_methylation_analysis/output/modified_bases.5mC.bed"
out_dir = "/home/ubuntu/ecoli_methylation_analysis/kmer_analysis"
possible_mods = pd.read_csv(megalodon,
                       names=["Chromosome", "start", "stop", "name", "score", "Strand", "start2", "end2", "color",
                              "coverage", "percentage"],
                       sep="\t",
                       usecols=["Chromosome", "start", "Strand", "coverage", "percentage"])

possible_mods['find'] = "C"
possible_mods['replace'] = "M"
possible_mods['ambig'] = "Y"

temp_mods = possible_mods[possible_mods["percentage"] != 100]
output_path = os.path.join(out_dir, f"megalodon_canonical.positions")
temp_mods.to_csv(output_path, sep="\t", index=False,
                 columns=["Chromosome", "start", "Strand", "find", "find"], header=False)

temp_mods = possible_mods[possible_mods["percentage"] != 0]
output_path = os.path.join(out_dir, f"megalodon_modified.positions")
temp_mods.to_csv(output_path, sep="\t", index=False,
                 columns=["Chromosome", "start", "Strand", "find", "replace"], header=False)

```

* Kmer breakdown
```
/home/ubuntu/functional_model_analysis/src/kmer_breakdown_of_positions_file.py --reference /home/ubuntu/ecoli_methylation_analysis/reference/ecoli.fa --positions_file /home/ubuntu/ecoli_methylation_analysis/kmer_analysis/megalodon.positions --output_dir /home/ubuntu/ecoli_methylation_analysis/kmer_analysis/

Number of M kmers: 6144
Number of found M kmers: 6136
Number of Canonical kmers: 3367
Number of found Canonical kmers: 3367
```

* Min site selection
```
import os
import pandas as pd
import numpy as np
from argparse import ArgumentParser
from py3helpers.seq_tools import ReferenceHandler, ReverseComplement
from py3helpers.utils import time_it, all_string_permutations

reference = "/home/ubuntu/ecoli_methylation_analysis/reference/ecoli.fa"
positions_file = "/home/ubuntu/ecoli_methylation_analysis/kmer_analysis/megalodon.positions"
output_dir = "/home/ubuntu/ecoli_methylation_analysis/kmer_analysis/"


positions_data = pd.read_csv(positions_file, names=["chr", "start", "strand", "find", "replace"], sep="\t")
positions_data["kmer"] = np.nan
# reference handler and reverse complement handler
rh = ReferenceHandler(reference)
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

# chr           U00096.3
# start             7481
# strand               +
# find                 C
# replace              M
# kmer       AACTCMGCTGA

sorted_dict = sorted(found_m_only_kmers.items(), key=lambda kv:(kv[1], kv[0]))
kmer_to_pos = {k: [] for k in all_string_permutations(alphabet, kmer_length)}
for i, x in mod_pos_data.iterrows():
    for i in range(kmer_length):
        k = x.kmer[i:kmer_length + i]
        if len(k) == kmer_length:
            kmer_to_pos[k] += [x]

top_kmer_pos = []
index = 0
for kmer, count in sorted_dict:
    pos_list = kmer_to_pos[kmer]
    if len(pos_list) == 0: 
        continue
    top_kmer_pos.append(pos_list[np.random.randint(len(pos_list))])

df = pd.DataFrame(top_kmer_pos)
# output_path = os.path.join(output_dir, f"top_megalodon_modified.positions")
# df.to_csv(output_path, sep="\t", index=False, columns=["chr", "start", "strand", "find", "replace"], header=False)

## canoncial
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

sorted_dict = sorted(found_no_m_kmers.items(), key=lambda kv:(kv[1], kv[0]))
kmer_to_pos = {k: [] for k in all_string_permutations(alphabet, kmer_length)}
for i, x in c_pos_data.iterrows():
    for i in range(kmer_length):
        k = x.kmer[i:kmer_length + i]
        if len(k) == kmer_length:
            kmer_to_pos[k] += [x]

# top_kmer_pos = []
index = 0
for kmer, count in sorted_dict:
    pos_list = kmer_to_pos[kmer]
    if len(pos_list) == 0: 
        continue
    top_kmer_pos.append(pos_list[np.random.randint(len(pos_list))])

df = pd.DataFrame(top_kmer_pos).drop_duplicates(subset=['chr', "start", "strand"])
output_path = os.path.join(output_dir, f"top_megalodon.positions")
df.to_csv(output_path, sep="\t", index=False, columns=["chr", "start", "strand", "find", "replace"], header=False)
```

```
import pandas as pd
import pyranges as pr
import os

delta = 6
output_dir = "/home/ubuntu/ecoli_methylation_analysis/kmer_analysis/"

output_path = os.path.join(output_dir, f"top_megalodon.positions")
df = pd.read_csv(output_path, names=["chr", "start", "strand", "find", "replace"], sep="\t")

df['Chromosome'] = df["chr"]
df['Strand'] = df["strand"]
df['Start'] = df['start'] - delta
df['End'] = df['start'] + delta

pyr = pr.PyRanges(df)
pyr = pyr.merge(strand=True)
pyr = pyr[(pyr.End - pyr.Start) == (delta * 2)]
result = pd.merge(pyr.df, df, how="left", on=["Chromosome", "Strand", "Start"], validate="one_to_one")
output_path = os.path.join(output_dir, f"final_top_megalodon.positions")
result.to_csv(output_path, sep="\t", index=False, columns=["chr", "start", "strand", "find", "replace"], header=False)
```

* Get reads
```
import pandas as pd

megalodon_txt = "/home/ubuntu/ecoli_methylation_analysis/output/per_read_modified_base_calls.txt"
positions_file = "/home/ubuntu/ecoli_methylation_analysis/kmer_analysis/final_top_megalodon.positions"
output_ids = "/home/ubuntu/ecoli_methylation_analysis/kmer_analysis/final_top_megalodon.readids.txt"

positions_data = pd.read_csv(positions_file, names=["chr", "start", "strand", "find", "replace"], sep="\t")


read_ids = []
targets = set([x[0] for x in positions_data.groupby(["start", "strand", "replace"])])
sub_targets = set([x[0] for x in positions_data.groupby(["start", "strand"])])
with open(megalodon_txt) as fh:
    fh.readline()
    for line in fh:
        if len(targets) == 0:
            break
        split_line = line.rstrip().split()
        read_id, chrm, strand, pos, mod_log_prob, can_log_prob, mod_base = split_line
        pos_strand = (int(pos), strand)
        if pos_strand in sub_targets:
            if float(mod_log_prob) < float(can_log_prob):
                call = "M"
            else:
                call = "C"
            pos_strand_call = (int(pos), strand, call)
            if pos_strand_call in targets:
                read_ids.append(read_id)
                targets -= {pos_strand_call}
                sub_targets -= {pos_strand}
                print(len(targets))
                
final_ids = set(read_ids)
with open(output_ids, "w") as fh:
    [print(x, file=fh) for x in final_ids]
                
```

