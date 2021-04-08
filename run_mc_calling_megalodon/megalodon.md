# Methylcalling using Megalodon
In order to get a better gauge on our accuracy vs the state of the art neural network based methylation calling 
algorithms, we decided to run Megalodon on our test reads.

## Source
https://github.com/nanoporetech/megalodon

### Procedure
megalodon raw_fast5s/ \
--outputs basecalls mappings mod_mappings mods \
--reference reference.fa --mod-motif m CG 0 \
--devices 0 1 --processes 40

