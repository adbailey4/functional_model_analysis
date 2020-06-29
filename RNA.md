# RNA analysis
* https://github.com/nanopore-wgs-consortium/NA12878/blob/master/nanopore-human-transcriptome/fastq_fast5_bulk.md

* Transcriptome
    * https://www.gencodegenes.org/human/release_27.html
    * ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.transcripts.fa.gz
    
* UBC run5
    * s3://nanopore-human-wgs/rna/UBC_Run5_20171031_DirectRNA
    * s3://nanopore-human-wgs/rna/fastq/UBC_Run5_20171031_DirectRNA.pass.dedup.fastq
    
* OICR run5
    * s3://nanopore-human-wgs/rna/OICR_Run5_20171101_DirectRNA
    * s3://nanopore-human-wgs/rna/fastq/OICR_Run5_20171101_DirectRNA.pass.dedup.fastq

* Positions
    * python /home/ubuntu/functional_model_analysis/src/create_canonical_rna_positions.py

* Models
    * ls | grep run-1 | grep rna | xargs -I{} python /home/ubuntu/functional_model_analysis/src/convert_csv_to_signalalign.py --dir /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/models/{} --rna --base_model /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/ont_canonical_rna.model --output_dir /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/rna_sa_models/{} --num_threads 16

* UBC
    * minimap2 --MD -t 10 -ax map-ont /home/ubuntu/mount/download/RNA_rel2/reference/gencode.v27.transcripts.fa /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/UBC_Run5_20171031_DirectRNA.pass.dedup.fastq | samtools view -@ 10 -bS - | samtools sort -@ 10 - > /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/UBC_Run5_20171031_DirectRNA.pass.dedup.sorted.bam && samtools view -@ 10 -bSF 2308 /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/UBC_Run5_20171031_DirectRNA.pass.dedup.sorted.bam > /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/UBC_Run5_20171031_DirectRNA.pass.dedup.2308.sorted.bam && samtools index /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/UBC_Run5_20171031_DirectRNA.pass.dedup.2308.sorted.bam
    * embed_main index -d /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fast5  /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/UBC_Run5_20171031_DirectRNA.pass.dedup.fastq
    * python /home/ubuntu/functional_model_analysis/src/get_reads_covering_positions.py --bam /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/UBC_Run5_20171031_DirectRNA.pass.dedup.2308.sorted.bam --positions_file /home/ubuntu/mount/download/RNA_rel2/reference/transcript_canonical.positions --output_file /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/canonical_covered_UBC_Run5_20171031_DirectRNA.pass.dedup.2308.sorted.bam
    * samtools merge /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/canonical_covered_UBC_Run5_20171031_DirectRNA.pass.dedup.2308.sorted.bam /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/canonical_covered_bams/canonical_covered_UBC_Run5_20171031_DirectRNA.pass.dedup.2308.sorted.bam1 /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/canonical_covered_bams/canonical_covered_UBC_Run5_20171031_DirectRNA.pass.dedup.2308.sorted.bam2 /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/canonical_covered_bams/canonical_covered_UBC_Run5_20171031_DirectRNA.pass.dedup.2308.sorted.bam3
    * python /home/ubuntu/functional_model_analysis/src/minimum_kmer_covering_reads.py --threads 16 --positions_file /home/ubuntu/mount/download/RNA_rel2/reference/transcript_canonical.positions --bam /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/canonical_covered_UBC_Run5_20171031_DirectRNA.pass.dedup.2308.sorted.bam --alphabet ACGT --kmer_length 5 --output_dir /home/ubuntu/mount/download/RNA_rel2/UBC_run5/UBC_canonical_read_coverage_stats --reference /home/ubuntu/mount/download/RNA_rel2/reference/gencode.v27.transcripts.fa
    * grep -f /home/ubuntu/mount/download/RNA_rel2/UBC_run5/UBC_canonical_read_coverage_stats/all_reads_covering_kmers_with_threshold_10.txt /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/UBC_Run5_20171031_DirectRNA.pass.dedup.fastq.index.readdb | grep fast5 > /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/cov_10_UBC_Run5_20171031_DirectRNA.pass.dedup.fastq.index.readdb
    * awk '{print $2}' /home/ubuntu/mount/download/RNA_rel2/UBC_run5/fastq/tmp2_cov_10_UBC_Run5_20171031_DirectRNA.pass.dedup.fastq.index.readdb | xargs -I{} cp {} /home/ubuntu/mount/download/RNA_rel2/UBC_run5/cov_10
    * sed 's='fast5/.*/fast5/.*/'='cov_10/'=g' 
    * python /home/ubuntu/functional_model_analysis/src/re_run_signalalign.py --dir /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/rna_sa_models/all_rna_models --output_dir /home/ubuntu/mount/UBC_runs/all_runs --base_model /home/ubuntu/mount/download/RNA_rel2/UBC_canonical_rna_runSignalAlign.config.json -v X --rna
    * find . -maxdepth 6 -name tempFiles_alignment | xargs tar --remove-files -czvf UBC_run_alignment_files.tar.gz
    * find . -maxdepth 6 -name variant_calls | xargs tar -czvf UBC_run_variant_calls.tar.gz

* OICR
    * minimap2 --MD -t 10 -ax map-ont /home/ubuntu/mount/download/RNA_rel2/reference/gencode.v27.transcripts.fa /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/OICR_Run5_20171101_DirectRNA.pass.dedup.fastq | samtools view -@ 10 -bS - | samtools sort -@ 10 - > /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/OICR_Run5_20171101_DirectRNA.pass.dedup.sorted.bam && samtools view -@ 10 -bSF 2308 /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/OICR_Run5_20171101_DirectRNA.pass.dedup.sorted.bam > /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/OICR_Run5_20171101_DirectRNA.pass.dedup.2308.sorted.bam && samtools index /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/OICR_Run5_20171101_DirectRNA.pass.dedup.2308.sorted.bam
    * embed_main index -d /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fast5  /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/OICR_Run5_20171101_DirectRNA.pass.dedup.fastq
    * python /home/ubuntu/functional_model_analysis/src/get_reads_covering_positions.py --bam /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/OICR_Run5_20171101_DirectRNA.pass.dedup.2308.sorted.bam --positions_file /home/ubuntu/mount/download/RNA_rel2/reference/transcript_canonical.positions --output_file /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/canonical_covered_bams/canonical_covered_OICR_Run5_20171101_DirectRNA.pass.dedup.2308.sorted.bam
    * samtools merge /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/canonical_covered_OICR_Run5_20171101_DirectRNA.pass.dedup.2308.sorted.bam /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/canonical_covered_bams/canonical_covered_OICR_Run5_20171101_DirectRNA.pass.dedup.2308.sorted.bam1 /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/canonical_covered_bams/canonical_covered_OICR_Run5_20171101_DirectRNA.pass.dedup.2308.sorted.bam2 /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/canonical_covered_bams/canonical_covered_OICR_Run5_20171101_DirectRNA.pass.dedup.2308.sorted.bam3
    * python /home/ubuntu/functional_model_analysis/src/minimum_kmer_covering_reads.py --threads 16 --positions_file /home/ubuntu/mount/download/RNA_rel2/reference/transcript_canonical.positions --bam /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/canonical_covered_OICR_Run5_20171101_DirectRNA.pass.dedup.2308.sorted.bam --alphabet ACGT --kmer_length 5 --output_dir /home/ubuntu/mount/download/RNA_rel2/OICR_run5/OICR_canonical_read_coverage_stats --reference /home/ubuntu/mount/download/RNA_rel2/reference/gencode.v27.transcripts.fa
    * grep -f /home/ubuntu/mount/download/RNA_rel2/OICR_run5/OICR_canonical_read_coverage_stats/all_reads_covering_kmers_with_threshold_10.txt /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/OICR_Run5_20171101_DirectRNA.pass.dedup.fastq.index.readdb | grep fast5 > /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/cov_10_OICR_Run5_20171101_DirectRNA.pass.dedup.fastq.index.readdb
    * awk '{print $2}' /home/ubuntu/mount/download/RNA_rel2/OICR_run5/fastq/tmp_cov_10_OICR_Run5_20171101_DirectRNA.pass.dedup.fastq.index.readdb | xargs -I{} cp {} /home/ubuntu/mount/download/RNA_rel2/OICR_run5/cov_10_fast5
    * sed 's='fast5/.*/fast5/.*/'='cov_10_fast5/'=g' 
    * python /home/ubuntu/functional_model_analysis/src/re_run_signalalign.py --dir /home/ubuntu/bisulfite_methylation_analysis/all_model_testing/rna_sa_models/all_rna_models --output_dir /home/ubuntu/mount/OICR_runs/all_runs --base_model /home/ubuntu/mount/download/RNA_rel2/OICR_canonical_rna_runSignalAlign.config.json -v X --rna
    * find . -maxdepth 6 -name tempFiles_alignment | xargs tar --remove-files -czvf OICR_run_alignment_files.tar.gz
    * find . -maxdepth 6 -name variant_calls | xargs tar -czvf OICR_run_variant_calls.tar.gz

