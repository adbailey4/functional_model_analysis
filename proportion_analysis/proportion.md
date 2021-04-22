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
```

* prep and upload SA data 
  ```bash /home/ubuntu/functional_model_analysis/proportion_analysis/prep_data.sh```
  
* Run signalalign using `run_multiple_sa.sh` and `proportion_job.yml`. See Makefile: prop.

