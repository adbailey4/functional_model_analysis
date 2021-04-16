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
  tar -xzf FAB39088/FAB39088_mc_not_cpg_fast5.tar.gz -C FAB39088/ && rm FAB39088/FAB39088_mc_not_cpg_fast5.tar.gz
  echo "FAB39088 Done"
  tar -xzf FAF01169/FAF01169_mc_not_cpg_fast5.tar.gz -C FAF01169/ && rm FAF01169/FAF01169_mc_not_cpg_fast5.tar.gz
  echo "FAF01169 Done"
  tar -xzf methyl_calling_all_C_position.tar.gz && rm methyl_calling_all_C_position.tar.gz
  echo "methyl_calling_all_C_position Done"
  sed -i 's='/data/'='$cwd/'=g' "$cwd"/FAB39088/mc_not_cpg_FAB39088.readdb
  sed -i 's='/data/'='$cwd/'=g' "$cwd"/FAF01169/mc_not_cpg_FAF01169.readdb
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
      "positions_file": null,
      "fast5_dirs": ["$cwd/FAB39088/mc_not_cpg_fast5"],
      "bwa_reference": "$cwd/all_C_methyl_replaced_references/forward.na12878_chr1.GRCh38_full_analysis_set_plus_decoy_hla.fa",
      "fofns": [],
      "readdb": "$cwd/FAB39088/mc_not_cpg_FAB39088.readdb",
      "fw_reference": "$cwd/all_C_methyl_replaced_references/forward.na12878_chr1.GRCh38_full_analysis_set_plus_decoy_hla.fa",
      "bw_reference": "$cwd/all_C_methyl_replaced_references/backward.na12878_chr1.GRCh38_full_analysis_set_plus_decoy_hla.fa",
      "kmers_from_reference": false,
      "motifs": null,
      "name": "FAB39088_na12878",
      "probability_threshold": 0.7,
      "number_of_kmer_assignments": 10000,
      "alignment_file": "$cwd/FAB39088/mc_not_cpg_FAB39088.2308.sorted.bam",
      "recursive": false,
      "assignments_dir": null
    },
    {
      "positions_file": null,
      "fast5_dirs": ["$cwd/FAF01169/mc_not_cpg_fast5_repack"],
      "bwa_reference": "$cwd/all_C_methyl_replaced_references/forward.na12878_chr1.GRCh38_full_analysis_set_plus_decoy_hla.fa",
      "fofns": [],
      "readdb": "$cwd/FAF01169/mc_not_cpg_FAF01169.readdb",
      "fw_reference": "$cwd/all_C_methyl_replaced_references/forward.na12878_chr1.GRCh38_full_analysis_set_plus_decoy_hla.fa",
      "bw_reference": "$cwd/all_C_methyl_replaced_references/backward.na12878_chr1.GRCh38_full_analysis_set_plus_decoy_hla.fa",
      "kmers_from_reference": false,
      "motifs": null,
      "name": "FAF01169_na12878",
      "probability_threshold": 0.7,
      "number_of_kmer_assignments": 10000,
      "alignment_file": "$cwd/FAF01169/mc_not_cpg_FAF01169.2308.sorted.bam",
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
embed_main sa2bed -d output/tempFiles_alignment/FAB39088_na12878/ -a ambig_model.model -o output/variant_calls/FAB39088_"$1".bed -t "$2" -c Y --overwrite
embed_main sa2bed -d output/tempFiles_alignment/FAF01169_na12878/ -a ambig_model.model -o output/variant_calls/FAF01169_"$1".bed -t "$2" -c Y --overwrite
rm config.json
}

main "$1" "$2" "$3"

