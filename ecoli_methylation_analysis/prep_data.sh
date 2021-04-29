n=20
# Select reads covering those positions
cat /home/ubuntu/ecoli_methylation_analysis/subset_reference/filtered_positions/top_"$n"_covered_filtered_* > /home/ubuntu/ecoli_methylation_analysis/subset_reference/filtered_positions/top_"$n"_variant.positions
python /home/ubuntu/functional_model_analysis/src/get_reads_covering_positions.py --bam /home/ubuntu/ecoli_methylation_analysis/fastq/MinION_EC_NAT.2308.sorted.bam --positions_file /home/ubuntu/ecoli_methylation_analysis/subset_reference/filtered_positions/top_"$n"_variant.positions --output_file /home/ubuntu/ecoli_methylation_analysis/sa_workspace/top_"$n"_variant.bam && mv /home/ubuntu/ecoli_methylation_analysis/sa_workspace/top_"$n"_variant.bam0 /home/ubuntu/ecoli_methylation_analysis/sa_workspace/top_"$n"_variant.bam

# Prep Data for SignalAlign
samtools view /home/ubuntu/ecoli_methylation_analysis/sa_workspace/top_"$n"_variant.bam | awk '{print $1}' > /home/ubuntu/ecoli_methylation_analysis/sa_workspace/top_"$n"_variant.readids.txt
grep -f /home/ubuntu/ecoli_methylation_analysis/sa_workspace/top_"$n"_variant.readids.txt /home/ubuntu/ecoli_methylation_analysis/fastq/MinION_EC_NAT.fastq.index.readdb | grep fast5 > /home/ubuntu/ecoli_methylation_analysis/sa_workspace/MinION_EC_NAT.fastq.index.readdb
rm /home/ubuntu/ecoli_methylation_analysis/sa_workspace/fast5/*
awk '{print $2}' /home/ubuntu/ecoli_methylation_analysis/sa_workspace/MinION_EC_NAT.fastq.index.readdb | xargs -I{} cp {} /home/ubuntu/ecoli_methylation_analysis/sa_workspace/fast5
sed -i 's='/home/ubuntu/ecoli_methylation_analysis/MinION_EC_NAT/.*/'='/home/ubuntu/ecoli_methylation_analysis/sa_workspace/fast5/'=g' /home/ubuntu/ecoli_methylation_analysis/sa_workspace/MinION_EC_NAT.fastq.index.readdb

# rename for consistency of downstream processes
mv /home/ubuntu/ecoli_methylation_analysis/sa_workspace/top_"$n"_variant.bam /home/ubuntu/ecoli_methylation_analysis/sa_workspace/ecoli_variant.bam
samtools sort /home/ubuntu/ecoli_methylation_analysis/sa_workspace/ecoli_variant.bam -o /home/ubuntu/ecoli_methylation_analysis/sa_workspace/ecoli_variant.sorted.bam
samtools index /home/ubuntu/ecoli_methylation_analysis/sa_workspace/ecoli_variant.sorted.bam
mv /home/ubuntu/ecoli_methylation_analysis/subset_reference/filtered_positions/top_"$n"_variant.positions /home/ubuntu/ecoli_methylation_analysis/sa_workspace/ecoli_variant.positions

# move to S3
tar -czf /home/ubuntu/ecoli_methylation_analysis/sa_workspace.tar.gz -C /home/ubuntu/ecoli_methylation_analysis/ sa_workspace/
aws s3 cp /home/ubuntu/ecoli_methylation_analysis/sa_workspace.tar.gz s3://bailey-k8s/functional_model_analysis/ecoli_methylation/ --profile bailey
