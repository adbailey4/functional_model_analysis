# functional_model_analysis
Information and scripts for functional analysis of kmer models


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