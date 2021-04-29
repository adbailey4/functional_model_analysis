#!/bin/bash
# $1 model bucket
# $2 job_count
# $3 output bucket

main() {
  # $1 model bucket
  # $2 job_count
  # $3 output bucket

  download_files
  run_sa_multiple_models "$1" "$2" "$3"
}

download_files() {
  cwd=$(pwd)
  echo "Extract Files"
  tar -xzf sa_workspace.tar.gz && rm sa_workspace.tar.gz
  sed -i 's='/home/ubuntu/ecoli_methylation_analysis/'='$cwd/'=g' "$cwd"/sa_workspace/MinION_EC_NAT.fastq.index.readdb
  echo "Download Done"
}

run_sa_multiple_models() {
  # $1 model bucket
  # $2 job_count
  # $3 output bucket
  mkdir models
  mkdir output
  for model in $(aws s3 ls s3://"$1")
  do
    if [[ "$model" == *.model ]]
    then
      aws s3 cp --no-progress s3://"$1""$model" models
      run_sa "$model" "$2"
      #               tar together
      tar -czf "$model".variant_calls.tar.gz -C output/ variant_calls/
      #              test upload
      aws s3 mv "$model".variant_calls.tar.gz s3://"$3"
      rm -r output/*
    fi
  done

}


run_sa() {
#             create output directory
# $1 model name
# $2 job_count
cwd=$(pwd)
MODEL_PATH="\"$cwd/models/$1\""

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
      "positions_file": "$cwd/sa_workspace/ecoli_variant.positions",
      "fast5_dirs": ["$cwd/sa_workspace/fast5"],
      "bwa_reference": "$cwd/ecoli.fa",
      "fofns": [],
      "readdb": "$cwd/sa_workspace/MinION_EC_NAT.fastq.index.readdb",
      "fw_reference": null,
      "bw_reference": null,
      "kmers_from_reference": false,
      "motifs": null,
      "name": "ecoli",
      "probability_threshold": 0.7,
      "number_of_kmer_assignments": 10000,
      "alignment_file": "$cwd/sa_workspace/ecoli_variant.sorted.bam",
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
  "output_dir": "output",
  "constraint_trim": null,
  "diagonal_expansion": null,
  "traceBackDiagonals": 150,
  "filter_reads": 0,
  "perform_kmer_event_alignment": true,
  "overwrite": true,
  "rna": false,
  "ambig_model": "$cwd/ambig_model.model"
}
EOF
#               run signalAlign and variant caller
echo "Running SignalAlign"
runSignalAlign.py run --config config.json
mkdir output/variant_calls
echo "Running sa2bed"
embed_main sa2bed -d output/tempFiles_alignment/ecoli/ -a ambig_model.model -o output/variant_calls/ecoli_"$1".bed -t "$2" -c Y --overwrite
rm config.json
}

main "$1" "$2" "$3"

