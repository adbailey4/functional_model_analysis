## Prepare Positions for Canonical and Methylation Calling


The NA12878 genome is commonly used as a benchmark for several DNA and RNA sequencing experiments. 
* ENCODE - "An integrated encyclopedia of DNA elements in the human genome). 
* Nanopore assembly and Methylation calling - "Nanopore sequencing and assembly of a human genome with ultra-long reads". 
<br />
<br />

Use bisulfite sequencing data from ENCODE (https://www.encodeproject.org/experiments/ENCSR890UQO/) to generate high confidence methylated cytosines and high confidence canonical bases.
<br />
<br />
<br />
<br />

### Process Bisulfite Data
#### Step1: Download

All data can be found https://www.encodeproject.org/experiments/ENCSR890UQO/.
* Two isogenic replicates.
    * Replicate 1: https://www.encodeproject.org/biosamples/ENCBS232JAP/
        * methylation state at CHG: 2.51 GB
            * https://www.encodeproject.org/files/ENCFF721BJM/@@download/ENCFF721BJM.bed.gz
        * methylation state at CpG: 636 MB
            * https://www.encodeproject.org/files/ENCFF279HCL/@@download/ENCFF279HCL.bed.gz
        * methylation at CHH: 7.98 GB
            * https://www.encodeproject.org/files/ENCFF448RTC/@@download/ENCFF448RTC.bed.gz
    * Replicate 2: https://www.encodeproject.org/biosamples/ENCBS232JAP/  
        * methylation state at CHG: 2.47 GB
            * https://www.encodeproject.org/files/ENCFF349NNL/@@download/ENCFF349NNL.bed.gz
        * methylation state at CpG: 635 MB
            * https://www.encodeproject.org/files/ENCFF835NTC/@@download/ENCFF835NTC.bed.gz
        * methylation at CHH: 7.88 GB
            * https://www.encodeproject.org/files/ENCFF038HXQ/@@download/ENCFF038HXQ.bed.gz
<br />
* Bed file format
    1. Reference chromosome or scaffold
    * Start position in chromosome
    * End position in chromosome
    * Name of item
    * Score from 0-1000. Capped number of reads
    * Strandedness, plus (+), minus (-), or unknown (.)
    * Start of where display should be thick (start codon)
    * End of where display should be thick (stop codon)
    * Color value (RGB)
    * Coverage, or number of reads
    * Percentage of reads that show methylation at this position in the genome


<br />

#### Step2: Compare and Aggregate
* Plot the correlation between these bed files in order to get an idea of our confidence levels of the "gold standard" of bisulfite sequencing. 
* Aggregate data
    * Create file with both CpG files aggregated together
    * Create file with both CHG files aggregated together
    * Create file with both CHH files aggregated together
    * Create file with CHH and CHG files aggregated together
    * Create file with all files aggregated together
    
#### Step3: Generate desired files
* Isolated Canonical Positions on Chr21
* CpG high confidence calls on Chr21
* mC high confidence calls on Chr21

    
    
<br />
<br />
<br />
Note: "Although DNA methylation occurs mostly in the CG context, it may also occur at CHG and CHH sites (where H can be any nucleotide other than G)." - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3903047/

