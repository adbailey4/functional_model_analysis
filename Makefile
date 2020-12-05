mmc:
	aws s3 cp mc_calling_multiple_models/run_multiple_sa.sh s3://bailey-k8s/functional_model_analysis/methyl_calling_data/
	# Run a kubernetes job with our image, prefix with USERNAME and timestamp
	TS=`date +"%Y%m%d-%H%M%S"` envsubst < mc_calling_multiple_models/mc_multiple_models_job.yml | kubectl create -f -

rdc:
	aws s3 cp run_dna_canonical/run_dna_canonical.sh s3://bailey-k8s/functional_model_analysis/canonical_dna_calling_data/
	# Run a kubernetes job with our image, prefix with USERNAME and timestamp
	TS=`date +"%Y%m%d-%H%M%S"` envsubst < run_dna_canonical/run_dna_canonical.yml | kubectl create -f -

