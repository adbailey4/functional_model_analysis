### Proportion analysis
Given the bisulfite sequencing dataset, we can select positions that are 
confidently modified at various fractions (0%, 10%, 20%, 30% .. 100%).

1) Select n positions of interest for each bin of fraction modified using bisulfite data.
2) Select m nanopore reads covering those positions.
3) Make predictions on nanopore reads

wget https://www.encodeproject.org/files/ENCFF038HXQ/@@download/ENCFF038HXQ.bed.gz && wget https://www.encodeproject.org/files/ENCFF448RTC/@@download/ENCFF448RTC.bed.gz