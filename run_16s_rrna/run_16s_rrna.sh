#!/bin/bash
# $1 model bucket
# $2 job_count
# $3 output bucket

main() {
  # $1 model bucket
  # $2 job_count
  # $3 output bucket
  DEBIAN_FRONTEND=noninteractive apt-get update
  DEBIAN_FRONTEND=noninteractive apt-get install -y awscli
  cp -r .aws ~/.aws
  #             download data
  run_sa_multiple_models "$1" "$2" "$3"
}


run_sa_multiple_models() {
  # $1 model bucket
  # $2 job_count
  # $3 output bucket
  mkdir /data/sa_analysis/models
  mkdir /data/sa_analysis/output
  for model in $(aws s3 ls s3://"$1" --profile bailey)
  do
    if [[ "$model" == *.model ]]
    then
      aws s3 cp --no-progress s3://"$1""$model" /data/sa_analysis/models/ --profile bailey
      run_sa "$model" "$2"
      #               tar together
      tar -czf "$model".variant_calls.tar.gz -C /data/sa_analysis/output/ variant_calls/
      #              upload
      aws s3 mv "$model".variant_calls.tar.gz s3://"$3" --profile bailey
      rm -r /data/sa_analysis/output/*
    fi
  done

}


run_sa() {
#             create output directory
# $1 model name
# $2 job_count
cwd=$(pwd)
MODEL_PATH="\"/data/sa_analysis/models/$1\""

echo "$cwd"
cat << EOF >> config.json
{
  "signal_alignment_args": {
    "target_regions": null,
    "track_memory_usage": false,
    "threshold": 0.01,
    "event_table": null,
    "embed": false,
    "delete_tmp": true,
    "output_format": "full"
  },
  "samples": [
    {
      "positions_file": "/data/reference/16S_final.positions",
      "fast5_dirs": ["/data/fast5"],
      "bwa_reference": "/data/reference/J01859.1.fa",
      "fofns": [],
      "readdb": "/data/sa_analysis/run_1_16s_reads.readdb",
      "fw_reference": null,
      "bw_reference": null,
      "kmers_from_reference": false,
      "motifs": null,
      "name": "ecoli_16S",
      "probability_threshold": 0.7,
      "number_of_kmer_assignments": 10000,
      "alignment_file": "/data/fastq/ecoli_16s_rrna_native.sorted.bam",
      "recursive": false,
      "assignments_dir": null
    }
  ],
  "path_to_bin": "/root/signalAlign/bin",
  "complement_hdp_model": null,
  "template_hdp_model": null,
  "complement_hmm_model": null,
  "template_hmm_model": $MODEL_PATH,
  "job_count": $2,
  "debug": false,
  "two_d": false,
  "output_dir": "/data/sa_analysis/output",
  "constraint_trim": null,
  "diagonal_expansion": null,
  "traceBackDiagonals": 150,
  "filter_reads": 0,
  "perform_kmer_event_alignment": true,
  "overwrite": true,
  "rna": true,
  "ambig_model": "/data/reference/ambig_chars.model"
}
EOF
#               run signalAlign and variant caller
echo "Running SignalAlign"
runSignalAlign.py run --config config.json
mkdir /data/sa_analysis/output/variant_calls
echo "Running sa2bed"
embed_main sa2bed -d /data/sa_analysis/output/tempFiles_alignment/ecoli_16S/ -a /data/reference/ambig_chars.model -o /data/sa_analysis/output/variant_calls/ecoli_16S_"$1".bed -t "$2" -c EF --overwrite
rm config.json
}

main "$1" "$2" "$3"

