# DNA Analysis
* https://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md

* Reference:
    * http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
    
* NOTTS: FAB39088
    * http://s3.amazonaws.com/nanopore-human-wgs/rel6/FASTQTars/FAB39088-288418386_Multi.tar
    * http://s3.amazonaws.com/nanopore-human-wgs/rel6/MultiFast5Tars/FAB39088-288418386_Multi_Fast5.tar

* Bham: FAF01169
    http://s3.amazonaws.com/nanopore-human-wgs/rel6/FASTQTars/FAF01169-4245879798_Multi.tar
    http://s3.amazonaws.com/nanopore-human-wgs/rel6/MultiFast5Tars/FAF01169-4245879798_Multi_Fast5.tar  

* Bisulfite Data:
    * [bisulfite](bisulfite.md)

* Positions
    * python /home/ubuntu/functional_model_analysis/src/create_canonical_positions.py

* Models
    * methyl
        * python /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/convert_csv_to_signalalign.py --dir /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/models/dna_native_r9.4_native-on-all_model --base_model /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/nanopolish_ACGMT_cpg.model --output_dir /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/sa_models/dna_native_r9.4_native-on-all_model --num_threads 48
    * canonical
        * ls | grep run-1 | grep dna | xargs -I{} python /home/ubuntu/functional_model_analysis/src/convert_csv_to_signalalign.py --dir /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/models/{} --base_model /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/ont_canonical.model --output_dir /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/sa_models/{} --num_threads 16

* FAB39088

    * embed_main split_multi_fast5 --fast5_dir /home/ubuntu/mount/download/FAB39088/Notts/FAB39088-288418386_Multi --output_dir /home/ubuntu/mount/download/FAB39088/fast5/ -j 8
    * embed_main index -d /home/ubuntu/mount/download/FAB39088/fast5/  /home/ubuntu/mount/download/FAB39088/fastq/FAB39088.fastq
    * minimap2 --MD -t 18 -ax map-ont /home/ubuntu/bisulfite_methylation_analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/mount/download/FAB39088/fastq/FAB39088.fastq | samtools view -@ 10 -bS - | samtools sort -@ 10 - > /home/ubuntu/mount/download/FAB39088/fastq/FAB39088.sorted.bam && samtools view -@ 40 -bSF 2308 /home/ubuntu/mount/download/FAB39088/fastq/FAB39088.sorted.bam > /home/ubuntu/mount/download/FAB39088/fastq/FAB39088.2308.sorted.bam && samtools index /home/ubuntu/mount/download/FAB39088/fastq/FAB39088.2308.sorted.bam
    * python /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/get_reads_covering_positions.py --bam /home/ubuntu/mount/download/FAB39088/fastq/FAB39088.2308.sorted.bam --positions_file /home/ubuntu/bisulfite_methylation_analysis/positions/canonical_cpg.positions --output_file /home/ubuntu/mount/download/FAB39088/fastq/canonical_cpg_FAB39088.2308.sorted.bam
    * python /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/get_reads_covering_positions.py --bam /home/ubuntu/mount/download/FAB39088/fastq/FAB39088.2308.sorted.bam --positions_file /home/ubuntu/bisulfite_methylation_analysis/positions/mc_not_cpg.positions --output_file /home/ubuntu/mount/download/FAB39088/fastq/mc_not_cpg_FAB39088.2308.sorted.bam
    * samtools view /home/ubuntu/mount/download/FAB39088/fastq/canonical_cpg_FAB39088.2308.sorted.bam | awk '{print $1}' > /home/ubuntu/mount/download/FAB39088/fastq/canonical_cpg_readids.txt
    * grep -f /home/ubuntu/mount/download/FAB39088/fastq/canonical_cpg_readids.txt /home/ubuntu/mount/download/FAB39088/fastq/FAB39088.fastq.index.readdb | grep fast5 > /home/ubuntu/mount/download/FAB39088/fastq/canonical_cpg.readdb2
    * awk '{print $2}' /home/ubuntu/mount/download/FAB39088/fastq/canonical_cpg.readdb | xargs -I{} cp {} /home/ubuntu/mount/download/FAB39088/canonical_fast5
    * samtools view /home/ubuntu/mount/download/FAB39088/fastq/mc_not_cpg_FAB39088.2308.sorted.bam | awk '{print $1}' > /home/ubuntu/mount/download/FAB39088/fastq/mc_not_cpg_readids.txt
    * grep -f /home/ubuntu/mount/download/FAB39088/fastq/mc_not_cpg_readids.txt /home/ubuntu/mount/download/FAB39088/fastq/FAB39088.fastq.index.readdb | grep fast5 > /home/ubuntu/mount/download/FAB39088/fastq/mc_not_cpg.readdb
    * awk '{print $2}' /home/ubuntu/mount/download/FAB39088/fastq/mc_not_cpg.readdb | xargs -I{} cp {} /home/ubuntu/mount/download/FAB39088/mc_not_cpg_fast5
    * sed 's='fast5/.*/'='canonical_fast5/'=g' 
    * sed 's='fast5/.*/'='mc_not_cpg_fast5/'=g' 
    * python /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/re_run_signalalign.py --dir /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/subsample_mc_calling/sa_models --output_dir /home/ubuntu/mount/FAB39088_runs/cpg_calling/subsample --base_model /home/ubuntu/mount/FAB39088_runs/mc_dna_runSignalAlign.config.json -v Y
    * python /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/re_run_signalalign.py --dir /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/sa_models/run-1_models --output_dir /home/ubuntu/mount/FAB39088_runs/canonical_calling/all_runs --base_model /home/ubuntu/mount/FAB39088_runs/canonical_dna_runSignalAlign.config.json -v X
    * find . -maxdepth 6 -name tempFiles_alignment | xargs tar --remove-files -czvf FAB39088_canonical_alignment_files.tar.gz
    * find . -maxdepth 6 -name variant_calls | xargs tar -czvf FAB39088_canonical_variant_calls.tar.gz
    * find . -maxdepth 6 -name tempFiles_alignment | xargs tar -czvf FAB39088_cpg_alignment_files.tar.gz
    * find . -maxdepth 6 -name variant_calls | xargs tar -czvf FAB39088_cpg_variant_calls.tar.gz

* FAF01169

    * split_multi_fast5 --fast5_dir /home/ubuntu/mount/download/FAF01169/Bham/FAF01169-4245879798_Multi --output_dir /home/ubuntu/mount/download/FAF01169/Bham/fast5/ -j 48
    * embed_main index -d /home/ubuntu/mount/download/FAF01169/Bham/fast5  /home/ubuntu/mount/download/FAF01169/Bham/fastq/FAF01169.fastq
    * minimap2 --MD -t 18 -ax map-ont /home/ubuntu/bisulfite_methylation_analysis/ref/GRCh38_full_analysis_set_plus_decoy_hla.fa /home/ubuntu/mount/download/FAF01169/Bham/fastq/FAF01169.fastq | samtools view -@ 10 -bS - | samtools sort -@ 10 - > /home/ubuntu/mount/download/FAF01169/Bham/fastq/FAF01169.2308.sorted.bam && samtools view -@ 40 -bSF 2308 /home/ubuntu/mount/download/FAF01169/Bham/fastq/FAF01169.2308.sorted.bam > /home/ubuntu/mount/download/FAF01169/Bham/fastq/FAF01169.2308.sorted.bam && samtools index /home/ubuntu/mount/download/FAF01169/Bham/fastq/FAF01169.2308.sorted.bam
    * python /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/get_reads_covering_positions.py --bam /home/ubuntu/mount/download/FAF01169/Bham/fastq/FAF01169.2308.sorted.bam --positions_file /home/ubuntu/bisulfite_methylation_analysis/positions/canonical_cpg.positions --output_file /home/ubuntu/mount/download/FAF01169/Bham/fastq/canonical_cpg_FAF01169.2308.sorted.bam
    * python /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/get_reads_covering_positions.py --bam /home/ubuntu/mount/download/FAF01169/Bham/fastq/FAF01169.2308.sorted.bam --positions_file /home/ubuntu/bisulfite_methylation_analysis/positions/mc_not_cpg.positions --output_file /home/ubuntu/mount/download/FAF01169/Bham/fastq/mc_not_cpg_FAF01169.2308.sorted.bam
    * samtools view /home/ubuntu/mount/download/FAF01169/Bham/fastq/canonical_cpg_FAF01169.2308.sorted.bam | awk '{print $1}' > /home/ubuntu/mount/download/FAF01169/Bham/fastq/canonical_cpg_readids.txt
    * grep -f /home/ubuntu/mount/download/FAF01169/Bham/fastq/canonical_cpg_readids.txt /home/ubuntu/mount/download/FAF01169/Bham/fastq/FAF01169.fastq.index.readdb | grep fast5 > /home/ubuntu/mount/download/FAF01169/Bham/fastq/canonical_cpg.readdb
    * awk '{print $2}' /home/ubuntu/mount/download/FAF01169/Bham/fastq/canonical_cpg.readdb | xargs -I{} cp {} /home/ubuntu/mount/download/FAF01169/Bham/canonical_fast5
    * samtools view /home/ubuntu/mount/download/FAF01169/Bham/fastq/mc_not_cpg_FAF01169.2308.sorted.bam | awk '{print $1}' > /home/ubuntu/mount/download/FAF01169/Bham/fastq/mc_not_cpg_readids.txt
    * grep -f /home/ubuntu/mount/download/FAF01169/Bham/fastq/mc_not_cpg_readids.txt /home/ubuntu/mount/download/FAF01169/Bham/fastq/FAF01169.fastq.index.readdb | grep fast5 > /home/ubuntu/mount/download/FAF01169/Bham/fastq/mc_not_cpg.readdb
    * awk '{print $2}' /home/ubuntu/mount/download/FAF01169/Bham/fastq/mc_not_cpg.readdb | xargs -I{} cp {} /home/ubuntu/mount/download/FAF01169/Bham/mc_not_cpg_fast5
    * sed 's='fast5/.*/'='canonical_fast5/'=g' 
    * sed 's='fast5/.*/'='mc_not_cpg_fast5/'=g' 
    * python /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/re_run_signalalign.py --dir /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/sa_models/dna_native_r9.4_native-on-all_model --output_dir /home/ubuntu/mount/FAB39088_runs/cpg_calling/all_runs --base_model /home/ubuntu/mount/FAB39088_runs/mc_dna_runSignalAlign.config.json -v Y 
    * python /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/re_run_signalalign.py --dir /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/sa_models/run-1_models --output_dir /home/ubuntu/mount/FAF01169_runs/canonical_calling/all_runs/ --base_model /home/ubuntu/mount/FAF01169_runs/canonical_dna_runSignalAlign.config.json -v X
    * find . -maxdepth 6 -name tempFiles_alignment | xargs tar --remove-files -czvf FAB39088_canonical_alignment_files.tar.gz
    * find . -maxdepth 6 -name variant_calls | xargs tar -czvf FAB39088_canonical_variant_calls.tar.gz
    * find . -maxdepth 6 -name tempFiles_alignment | xargs tar -czvf FAB39088_cpg_alignment_files.tar.gz
    * find . -maxdepth 6 -name variant_calls | xargs tar -czvf FAB39088_cpg_variant_calls.tar.gz

