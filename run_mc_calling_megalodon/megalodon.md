# Methylcalling using Megalodon
In order to get a better gauge on our accuracy vs the state of the art neural network based methylation calling 
algorithms, we decided to run Megalodon on our test reads.

## Source
https://github.com/nanoporetech/megalodon

### Pipeline
We needed to re-generate the fast5 files because the in house fast5 
splitting tool we use is not as good as the ONT version (https://github.com/nanoporetech/ont_fast5_api).
So, we downloaded the raw Fast5s and split using multi_to_single_fast5. 
After we split the fast5s, we isolated the reads of interest using the following python code.

```
from py3helpers.utils import list_dir_recursive
import pandas as pd
import os
from shutil import copy

# FAB
new_split_dir = "/home/ubuntu/predictive_analysis/dna_megalodon/temp_data/Notts/FAB39088_split_fast5"
readdb = "/home/ubuntu/predictive_analysis/dna_megalodon/data/FAB39088/mc_not_cpg_FAB39088.readdb"
output = "/home/ubuntu/predictive_analysis/dna_megalodon/new_split/FAB39088/split_fast5"
# FAF
new_split_dir = "/home/ubuntu/predictive_analysis/dna_megalodon/temp_data/Bham/FAF01169_split_fast5"
readdb = "/home/ubuntu/predictive_analysis/dna_megalodon/data/FAF01169/mc_not_cpg_FAF01169.readdb"
output = "/home/ubuntu/predictive_analysis/dna_megalodon/new_split/FAF01169/split_fast5"


df_r = pd.read_csv(readdb, sep="\t", header=None)
read_ids = set(df_r[0])
files = []
for x in list_dir_recursive(new_split_dir, ext="fast5"):
    name = os.path.splitext(os.path.split(x)[1])[0]
    if name in read_ids:
        files.append(x)
        copy(x, output)
    
print(len(files))
print(len(read_ids))

```

### megalodon commands

megalodon /home/ubuntu/predictive_analysis/dna_megalodon/new_split/FAF01169/split_fast5 \
--outputs mod_mappings mods \
--reference /home/ubuntu/predictive_analysis/dna_megalodon/data/GRCh38_full_analysis_set_plus_decoy_hla.fa \
--processes 1 --overwrite --guppy-server-path guppy_basecall_server \
--output-directory /home/ubuntu/predictive_analysis/dna_megalodon/new_split/FAF01169/megalodon_output \
--guppy-timeout 1000 --guppy-concurrent-reads 1 --guppy-params "--num_callers 7 --cpu_threads_per_caller 10 --chunks_per_runner 100"

[comment]: <> (--chunks_per_runner 1)

megalodon /home/ubuntu/predictive_analysis/dna_megalodon/new_split/FAB39088/split_fast5 \
--outputs mod_mappings mods \
--reference /home/ubuntu/predictive_analysis/dna_megalodon/data/GRCh38_full_analysis_set_plus_decoy_hla.fa \
--processes 96 --overwrite --guppy-server-path guppy_basecall_server \
--output-directory /home/ubuntu/predictive_analysis/dna_megalodon/new_split/FAB39088/megalodon_output

megalodon_extras per_read_text modified_bases .

#### Local Calling
Local calling on my own GPU requires more CPU memory for loading the reference and decreasing the GPU memory footprint

megalodon /home/bailey/data/new_split/FAF01169/split_fast5 \
--outputs mod_mappings mods \
--reference /home/bailey/data/GRCh38_full_analysis_set_plus_decoy_hla.fa \
--devices cuda:0 --processes 1  --overwrite --guppy-server-path guppy_basecall_server \
--output-directory /home/bailey/data/new_split/FAF01169/megalodon_output \
--guppy-timeout 500 --guppy-concurrent-reads 1 --guppy-config dna_r9.4.1_450bps_modbases_5mc_hac.cfg \
--guppy-params "--gpu_runners_per_device 1 --chunks_per_runner 1"






### Convert text outputs to same format as our other callers

```
import pandas as pd
import numpy as np

text_output = "/home/ubuntu/predictive_analysis/dna_megalodon/new_split/FAB39088/megalodon_output/per_read_modified_base_calls.txt"
text_output = "/home/ubuntu/predictive_analysis/dna_megalodon/new_split/FAF01169/megalodon_output/per_read_modified_base_calls.txt"
positions = "/home/ubuntu/predictive_analysis/all_C.positions"

positions_data = pd.read_csv(positions, sep='\t', header=None)
mt_data = pd.read_csv(text_output, sep='\t')

positions_data['combined'] = positions_data[[0,1,2]].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)
mt_data['combined'] = mt_data[["chrm", "pos", "strand"]].apply(lambda row: '_'.join(row.values.astype(str)), axis=1)

mask = mt_data['combined'].isin(positions_data['combined'])

columns = ["read_id","contig","reference_index","strand","variants","prob1","prob2","prob3","prob4"]
#06de8374-698e-4bef-8501-781c0eecb41e,chr5,44931258,-,CM,0.93582,0.0641801,,
new_data = pd.DataFrame(columns=columns)
new_data["read_id"] = megalodon_text_data["read_id"]
new_data["contig"] = megalodon_text_data["chrm"]
new_data["reference_index"] = megalodon_text_data["pos"]
new_data["strand"] = megalodon_text_data["strand"]
new_data["variants"] = CM
new_data["prob1"] = megalodon_text_data["mod_log_prob"]
new_data["prob2"] = megalodon_text_data["can_log_prob"]

```
