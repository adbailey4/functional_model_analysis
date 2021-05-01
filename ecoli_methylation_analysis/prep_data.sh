# Prep Data for SignalAlign
samtools view /home/ubuntu/ecoli_methylation_analysis/fastq/MinION_EC_NAT.2308.sorted.bam | grep -f /home/ubuntu/ecoli_methylation_analysis/kmer_analysis/final_top_megalodon.readids.txt > /home/ubuntu/ecoli_methylation_analysis/sa_workspace/final_top_megalodon_noheader.sam
cat /home/ubuntu/ecoli_methylation_analysis/fastq/header.txt /home/ubuntu/ecoli_methylation_analysis/sa_workspace/final_top_megalodon_noheader.sam > /home/ubuntu/ecoli_methylation_analysis/sa_workspace/final_top_megalodon.sam
samtools view -Sb /home/ubuntu/ecoli_methylation_analysis/sa_workspace/final_top_megalodon.sam -o /home/ubuntu/ecoli_methylation_analysis/sa_workspace/final_top_megalodon.bam

grep -f /home/ubuntu/ecoli_methylation_analysis/kmer_analysis/final_top_megalodon.readids.txt /home/ubuntu/ecoli_methylation_analysis/fastq/MinION_EC_NAT.fastq.index.readdb | grep fast5 > /home/ubuntu/ecoli_methylation_analysis/sa_workspace/MinION_EC_NAT.fastq.index.readdb
rm /home/ubuntu/ecoli_methylation_analysis/sa_workspace/fast5/*
awk '{print $2}' /home/ubuntu/ecoli_methylation_analysis/sa_workspace/MinION_EC_NAT.fastq.index.readdb | xargs -I{} cp {} /home/ubuntu/ecoli_methylation_analysis/sa_workspace/fast5
sed -i 's='/home/ubuntu/ecoli_methylation_analysis/MinION_EC_NAT/.*/'='/home/ubuntu/ecoli_methylation_analysis/sa_workspace/fast5/'=g' /home/ubuntu/ecoli_methylation_analysis/sa_workspace/MinION_EC_NAT.fastq.index.readdb

# rename for consistency of downstream processes
mv /home/ubuntu/ecoli_methylation_analysis/sa_workspace/final_top_megalodon.bam /home/ubuntu/ecoli_methylation_analysis/sa_workspace/ecoli_variant.bam
samtools sort /home/ubuntu/ecoli_methylation_analysis/sa_workspace/ecoli_variant.bam -o /home/ubuntu/ecoli_methylation_analysis/sa_workspace/ecoli_variant.sorted.bam
samtools index /home/ubuntu/ecoli_methylation_analysis/sa_workspace/ecoli_variant.sorted.bam

cp /home/ubuntu/ecoli_methylation_analysis/kmer_analysis/final_top_megalodon.positions /home/ubuntu/ecoli_methylation_analysis/sa_workspace/ecoli_variant.positions
sed -i 's/C$/Y/g; s/M$/Y/g' ecoli_variant.positions

# move to S3
tar -czf /home/ubuntu/ecoli_methylation_analysis/sa_workspace.tar.gz -C /home/ubuntu/ecoli_methylation_analysis/ sa_workspace/
aws s3 cp /home/ubuntu/ecoli_methylation_analysis/sa_workspace.tar.gz s3://bailey-k8s/functional_model_analysis/ecoli_methylation/ --profile bailey
