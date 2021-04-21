### Proportion analysis
Given the bisulfite sequencing dataset, we can select positions that are
confidently modified at various fractions (10%, 20%, 30% .. 90%).

1) Select n positions of interest for each bin of fraction modified using bisulfite data.
2) Select m nanopore reads covering those positions.
3) Make predictions on nanopore reads

### Proportion analysis Workflow

* Filtering positions
  Within Chromosome 1, if if both rep1 and rep2 have >10 reads and are +-1% around 0, 10 ... 100% then we write out those positions.
```
grep -P 'chr1\t' ENCFF038HXQ.bed > chr1_ENCFF038HXQ.bed
grep -P 'chr1\t' ENCFF349NNL.bed > chr1_ENCFF349NNL.bed
grep -P 'chr1\t' ENCFF721BJM.bed > chr1_ENCFF721BJM.bed
grep -P 'chr1\t' ENCFF279HCL.bed > chr1_ENCFF279HCL.bed
grep -P 'chr1\t' ENCFF448RTC.bed > chr1_ENCFF448RTC.bed
grep -P 'chr1\t' ENCFF835NTC.bed > chr1_ENCFF835NTC.bed
python create_proportion_analysis_positions.py  
```

* Align all reads to chromosome 1
```
minimap2 --MD -t 40 -ax map-ont /home/ubuntu/bisulfite_methylation_analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/bisulfite_methylation_analysis/fastq/FAB39088.fastq | samtools view -@ 8 -bS - | samtools sort -@ 8 - > /home/ubuntu/bisulfite_methylation_analysis/fastq/FAB39088.sorted.bam && samtools view -@ 8 -bSF 2308 /home/ubuntu/bisulfite_methylation_analysis/fastq/FAB39088.sorted.bam > /home/ubuntu/bisulfite_methylation_analysis/fastq/FAB39088.2308.sorted.bam && samtools index /home/ubuntu/bisulfite_methylation_analysis/fastq/FAB39088.2308.sorted.bam
minimap2 --MD -t 40 -ax map-ont /home/ubuntu/bisulfite_methylation_analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/bisulfite_methylation_analysis/fastq/FAF01169.fastq | samtools view -@ 8 -bS - | samtools sort -@ 8 - > /home/ubuntu/bisulfite_methylation_analysis/fastq/FAF01169.sorted.bam && samtools view -@ 8 -bSF 2308 /home/ubuntu/bisulfite_methylation_analysis/fastq/FAF01169.sorted.bam > /home/ubuntu/bisulfite_methylation_analysis/fastq/FAF01169.2308.sorted.bam && samtools index /home/ubuntu/bisulfite_methylation_analysis/fastq/FAF01169.2308.sorted.bam

python /home/ubuntu/functional_model_analysis/src/get_reads_covering_positions.py --bam /home/ubuntu/bisulfite_methylation_analysis/fastq/FAF01169.2308.sorted.bam --positions_file /home/ubuntu/bisulfite_methylation_analysis/positions/filtered_chr1_all.positions --output_file /home/ubuntu/bisulfite_methylation_analysis/fastq/FAF01169_filtered_chr1_all.bam
python /home/ubuntu/functional_model_analysis/src/get_reads_covering_positions.py --bam /home/ubuntu/bisulfite_methylation_analysis/fastq/FAB39088.2308.sorted.bam --positions_file /home/ubuntu/bisulfite_methylation_analysis/positions/filtered_chr1_all.positions --output_file /home/ubuntu/bisulfite_methylation_analysis/fastq/FAB39088_filtered_chr1_all.bam
```

* Select top N most covered positions for each fraction

```
python /home/ubuntu/functional_model_analysis/proportion_analysis/top_n_most_covered_positions.py --bam /home/ubuntu/bisulfite_methylation_analysis/fastq/FAB39088_filtered_chr1_all.sorted.bam --positions_file /home/ubuntu/bisulfite_methylation_analysis/positions/filtered_chr1_10_variant.positions --output_dir /home/ubuntu/bisulfite_methylation_analysis/filtered_positions -n 100 --positions_directory /home/ubuntu/bisulfite_methylation_analysis/variant_positions
```
* Select reads covering those positions

```
cat top_100_covered_filtered_chr1_* > top_100_chr1_variant.positions
python /home/ubuntu/functional_model_analysis/src/get_reads_covering_positions.py --bam /home/ubuntu/bisulfite_methylation_analysis/fastq/FAF01169.2308.sorted.bam --positions_file /home/ubuntu/bisulfite_methylation_analysis/filtered_positions/top_100_chr1_variant.positions --output_file /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_100_chr1_variant.bam && mv /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_100_chr1_variant.bam0 /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_100_chr1_variant.bam
python /home/ubuntu/functional_model_analysis/src/get_reads_covering_positions.py --bam /home/ubuntu/bisulfite_methylation_analysis/fastq/FAB39088.2308.sorted.bam --positions_file /home/ubuntu/bisulfite_methylation_analysis/filtered_positions/top_100_chr1_variant.positions --output_file /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_100_chr1_variant.bam && mv /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_100_chr1_variant.bam0 /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_100_chr1_variant.bam
```

* Prep Data for SignalAlign

```
samtools view /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_100_chr1_variant.bam | awk '{print $1}' > /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_100_chr1_variant_readids.txt
samtools view /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_100_chr1_variant.bam | awk '{print $1}' > /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_100_chr1_variant_readids.txt
grep -f /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_100_chr1_variant_readids.txt /home/ubuntu/bisulfite_methylation_analysis/fastq/FAF01169.fastq.index.readdb | grep fast5 > /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169.readdb
grep -f /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_100_chr1_variant_readids.txt /home/ubuntu/bisulfite_methylation_analysis/fastq/FAB39088.fastq.index.readdb | grep fast5 > /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088.readdb
awk '{print $2}' /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088.readdb | xargs -I{} cp {} /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/fast5
awk '{print $2}' /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169.readdb | xargs -I{} cp {} /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/fast5
sed -i 's='/home/ubuntu/predictive_analysis/dna_megalodon/data/all_data/Notts/FAB39088_split_fast5/.*/'='/data/FAB39088/fast5/'=g' /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088.readdb
sed -i 's='/home/ubuntu/predictive_analysis/dna_megalodon/data/all_data/Bham/FAF01169_split_fast5/.*/'='/data/FAF01169/fast5/'=g' /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169.readdb
```



