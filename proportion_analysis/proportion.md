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

Align all reads to chromosome 1 and get coverage counts for each position.
Select top N most covered positions for each fraction
Select reads covering those positions
Run SignalAlign


