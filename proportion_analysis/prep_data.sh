n=200
python /home/ubuntu/functional_model_analysis/proportion_analysis/top_n_most_covered_positions.py --bam /home/ubuntu/bisulfite_methylation_analysis/fastq/FAF01169.2308.sorted.bam --output_dir /home/ubuntu/bisulfite_methylation_analysis/filtered_positions/FAF -n "$n" --positions_directory /home/ubuntu/bisulfite_methylation_analysis/variant_positions
python /home/ubuntu/functional_model_analysis/proportion_analysis/top_n_most_covered_positions.py --bam /home/ubuntu/bisulfite_methylation_analysis/fastq/FAB39088.2308.sorted.bam --output_dir /home/ubuntu/bisulfite_methylation_analysis/filtered_positions/FAB -n "$n" --positions_directory /home/ubuntu/bisulfite_methylation_analysis/variant_positions

# Select reads covering those positions
cat /home/ubuntu/bisulfite_methylation_analysis/filtered_positions/FAB/top_"$n"_covered_filtered_chr1_* > /home/ubuntu/bisulfite_methylation_analysis/filtered_positions/top_"$n"_chr1_variant.positions
cat /home/ubuntu/bisulfite_methylation_analysis/filtered_positions/FAF/top_"$n"_covered_filtered_chr1_* >> /home/ubuntu/bisulfite_methylation_analysis/filtered_positions/top_"$n"_chr1_variant.positions

python /home/ubuntu/functional_model_analysis/src/get_reads_covering_positions.py --bam /home/ubuntu/bisulfite_methylation_analysis/fastq/FAF01169.2308.sorted.bam --positions_file /home/ubuntu/bisulfite_methylation_analysis/filtered_positions/top_"$n"_chr1_variant.positions --output_file /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_"$n"_chr1_variant.bam && mv /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_"$n"_chr1_variant.bam0 /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_"$n"_chr1_variant.bam
python /home/ubuntu/functional_model_analysis/src/get_reads_covering_positions.py --bam /home/ubuntu/bisulfite_methylation_analysis/fastq/FAB39088.2308.sorted.bam --positions_file /home/ubuntu/bisulfite_methylation_analysis/filtered_positions/top_"$n"_chr1_variant.positions --output_file /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_"$n"_chr1_variant.bam && mv /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_"$n"_chr1_variant.bam0 /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_"$n"_chr1_variant.bam

# Prep Data for SignalAlign

samtools view /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_"$n"_chr1_variant.bam | awk '{print $1}' > /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_"$n"_chr1_variant_readids.txt
samtools view /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_"$n"_chr1_variant.bam | awk '{print $1}' > /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_"$n"_chr1_variant_readids.txt
grep -f /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_"$n"_chr1_variant_readids.txt /home/ubuntu/bisulfite_methylation_analysis/fastq/FAF01169.fastq.index.readdb | grep fast5 > /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169.readdb
grep -f /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_"$n"_chr1_variant_readids.txt /home/ubuntu/bisulfite_methylation_analysis/fastq/FAB39088.fastq.index.readdb | grep fast5 > /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088.readdb
rm /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/fast5/*
rm /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/fast5/*
awk '{print $2}' /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088.readdb | xargs -I{} cp {} /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/fast5
awk '{print $2}' /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169.readdb | xargs -I{} cp {} /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/fast5
sed -i 's='/home/ubuntu/predictive_analysis/dna_megalodon/data/all_data/Notts/FAB39088_split_fast5/.*/'='/data/FAB39088/fast5/'=g' /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088.readdb
sed -i 's='/home/ubuntu/predictive_analysis/dna_megalodon/data/all_data/Bham/FAF01169_split_fast5/.*/'='/data/FAF01169/fast5/'=g' /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169.readdb

# rename for consistency of downstream processes
mv /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_"$n"_chr1_variant.bam /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_chr1_variant.bam
samtools sort /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_chr1_variant.bam -o /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_chr1_variant.sorted.bam
samtools index /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169/FAF01169_top_chr1_variant.sorted.bam
mv /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_"$n"_chr1_variant.bam /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_chr1_variant.bam
samtools sort /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_chr1_variant.bam -o /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_chr1_variant.sorted.bam
samtools index /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088/FAB39088_top_chr1_variant.sorted.bam


mv /home/ubuntu/bisulfite_methylation_analysis/filtered_positions/top_"$n"_chr1_variant.positions /home/ubuntu/bisulfite_methylation_analysis/filtered_positions/top_chr1_variant.positions

tar -czf /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169.tar.gz -C /home/ubuntu/bisulfite_methylation_analysis/data/ FAF01169/
tar -czf /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088.tar.gz -C /home/ubuntu/bisulfite_methylation_analysis/data/ FAB39088/

aws s3 cp /home/ubuntu/bisulfite_methylation_analysis/data/FAF01169.tar.gz s3://bailey-k8s/functional_model_analysis/proportional_mc_data/ --profile bailey
aws s3 cp /home/ubuntu/bisulfite_methylation_analysis/data/FAB39088.tar.gz s3://bailey-k8s/functional_model_analysis/proportional_mc_data/ --profile bailey
aws s3 cp /home/ubuntu/bisulfite_methylation_analysis/filtered_positions/top_chr1_variant.positions s3://bailey-k8s/functional_model_analysis/proportional_mc_data/ --profile bailey
