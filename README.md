# functional_model_analysis
Information and scripts for functional analysis of kmer models

*  Required programs (in PATH)
    * signalAlign 
        * https://github.com/UCSC-nanopore-cgl/signalAlign: commit 3baaf3e47536278d636f74d2afadf30f1304ee4e
    * embed_fast5 
        * https://github.com/adbailey4/embed_fast5: commit 70c7f1ff009615c9c9e2612c49909cc74a966f55
    * samtools
    * minimap2
    * bwa
    
* Workflow
    * [DNA Analysis](DNA.md)
    * [RNA Analysis](RNA.md)
    * [ribosomal RNA analysis](rrna.md)

#### Docker
* The most accurate way to reproduce results is to use the docker image adbailey4/signalalign@sha256:a350ce89a00e23b96f2224a0ca8fc84e53ba7d44fde2c75331218c73b4833b1a for your analysis
* All final signalAlign runs were submitted through a kubernetes cluster using the job configs and bash scripts located in [mc_calling_multiple_models](mc_calling_multiple_models) and [run_dna_canonical](run_dna_canonical) 

## Scripts
#### convert_csv_to_signalalign.py
* This script was used to convert the output format of the machine learning model to the signalAlign model format

```text
Convert csv model to signalalign model

optional arguments:
  -h, --help            show this help message and exit
  --dir DIR, -d DIR     Path to directory of csv files to convert to
                        signalalign models
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Path to directory to output signalalign models
  --base_model BASE_MODEL, -m BASE_MODEL
                        Path to base signalalign model
  --num_threads NUM_THREADS, -t NUM_THREADS
                        Number of threads to run
  --rna, -r             Boolean option if models are for rna
```

#### create_canonical_positions.py
* This is the script used to read in all NA12878 bisulfite sequencing data and generate corresponding confident canonical and 5mC sites.  

#### create_canonical_rna_positions.py
* This is the script used to create canonical RNA positions.  

#### kmer_breakdown_of_positions_file.py
* Script to determine kmer counts and totals of a positions file of M and C positions

```
Write kmer breakdown by position given a positions file and reference

optional arguments:
  -h, --help            show this help message and exit
  --reference REFERENCE, -r REFERENCE
                        Path to reference
  --positions_file POSITIONS_FILE, -p POSITIONS_FILE
                        Path to positions_file
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Path to output directory

```

#### coverage_from_variant_calls.py
* Get coverage from RNA canonical, DNA canonical and DNA methylation variant call files to get final kmer coverage counts.  
* The following was used to calculate specific kmer breakdowns of the methyl kmer analysis
```python
file_path = "/home/ubuntu/mount/kmer_counts/FAB39088_methyl.tsv"
file_path = "/home/ubuntu/mount/kmer_counts/FAF01169_methyl.tsv"
file_path = "/home/ubuntu/mount/kmer_counts/FAB39088_methyl_FAF01169_methyl.tsv"

with open(file_path, 'r') as fh:
    lines = fh.readlines()
    data = [x.rstrip().split() for x in lines]
    m_counter = 0
    m_counts = 0
    m_coverage = 0 
    c_counter = 0
    c_counts = 0
    c_coverage = 0 
    
    for x in data:
        x_cov = float(x[1])
        if x[0].count("M") == 1:
            m_counter += 1
            if x_cov > 0:
                m_counts += 1
                m_coverage += x_cov
        if x[0].count("M") == 0 and x[0].count("C") >= 1:
            c_counter += 1
            if x_cov > 0:
                c_counts += 1
                c_coverage += x_cov
    
    print(m_counts, m_counter, m_coverage/m_counts)
    print(c_counts, c_counter, c_coverage/c_counts)
```
#### get_reads_covering_positions.py
* Simple script which selects reads which cover specific positions. This is used to narrow down the search for reads covering confident canonical and 5mC sites.

```text
Get all reads covering positions file

optional arguments:
  -h, --help            show this help message and exit
  --bam BAM, -d BAM     Bam file to look through
  --positions_file POSITIONS_FILE, -p POSITIONS_FILE
                        Path to base positions file
  --output_file OUTPUT_FILE, -o OUTPUT_FILE
                        Path to new bam file

```

#### mininum_kmer_covering_reads.py

* Get kmer coverage information and select for the fewest number of reads which cover the most number of kmers

```text
Get the minimum number of reads that cover all kmers given a positions file

optional arguments:
  -h, --help            show this help message and exit
  --bam BAM, -b BAM     Path to bam file
  --positions_file POSITIONS_FILE, -p POSITIONS_FILE
                        Path to positions_file
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Path to output directory
  --alphabet ALPHABET, -a ALPHABET
                        Set of characters to cover
  --kmer_length KMER_LENGTH, -l KMER_LENGTH
                        Length of the kmer to cover
  --reference REFERENCE, -r REFERENCE
                        Reference if MD flag is not set in BAM

```

#### re_run_signalalign.py

* run signalAlign multiple times for each model in a given directory

```text
Re-run signalalign multiple times for each model

optional arguments:
  -h, --help            show this help message and exit
  --dir DIR, -d DIR     Path to directory of model files
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Path to directory to output signalalign runs
  --base_model BASE_MODEL, -m BASE_MODEL
                        Base model file for signalalign
  --variants VARIANTS, -v VARIANTS
                        variants to analyze
  --rna                 set if rna reads

```


### Summary Of Event Data

Considering the computational complexity of signalAlign, for selected prediction sites as mentioned above, we performed the following filtering steps on NA12878 cell line native genomic DNA and mRNA datasets to get the most efficient read sets. First, we calculated read-level kmer coverage. For example, the center T-site of DNA read CAGATTACAGA was selected for signalAlign prediction. 6mers CAGATT, AGATTA, GATTAC, ATTACA, TTACAG and TACAGA span such T-site, therefore considered as being covered. Based on such read-level kmer coverage, we iteratively selected reads containing the least frequently covered kmers, building the most efficient read sets. For canonical DNA functional analysis, our final FAB39088 set contained 336 reads, which covered 3692/4096 canonical DNA 6mers with an average 13.62x coverage. And the final FAF01169 set contained 339 reads, which covered 3717/4096 canonical DNA 6mers with an average 18.10x coverage. Combining the two sets, the total average coverage was 31.72x for 3886/4096 kmers. For M DNA functional analysis, our final FAB39088 set contained 1706 reads, which covered 2625/3367 C-only DNA 6mers with an average 61.52x coverage as negative control, and 3105/6144 single-M DNA 6mers with an average 5.01x coverage. The final FAF01169 set contained 1396 reads, which covered 2610/3367 C-only DNA 6mers with an average 63.26x coverage as negative control, and 3140/6144 single-M DNA 6mers with an average 4.76x coverage. Combining the two sets, in total 2792/3367 C-only DNA 6mers were covered with an average  58.49x coverage, and 3481/6144 single-M DNA 6mers were covered with an average 4.38x coverage. FAB39088 and FAF01169 denoted the two biological replicates of NA12878 cell line native genomic DNA sequencing experiments. For canonical RNA functional analysis, our final UBC set contained 146 reads, which covered 1024/1024 canonical RNA 5mers with an average 11.59 coverage. And the final OICR set contained 163 reads, which covered 1024/1024 canonical RNA 5mers with an average 13.13 coverage. Combining the two sets, in total 1024/1024 canonical RNA 5mers were covered with an average 24.72 coverage. UBC and OICR denote the two biological replicates of NA12878 cell line native mRNA sequencing experiments.
