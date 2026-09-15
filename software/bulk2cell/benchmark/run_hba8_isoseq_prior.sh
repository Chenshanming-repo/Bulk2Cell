#!/usr/bin/env bash
# Full HBA8: Iso-Seq abundance prior only; append --resume for this output.
set -euo pipefail
cd /homeb/user/AD_data/AD-MultiOmics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
exec taskset --cpu-list 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29 bulk2cell/.venv/bin/python bulk2cell/benchmark/run_full_fit_workers.py --fit-workers 20 -- quantify-full --catalog /homeb/user/AD_data/AD-MultiOmics/bulk2cell/results/HBA8_full/inputs/catalog.json.gz --bam /homeb/user/AD_data/AD-MultiOmics/bulk2cell/results/HBA8_full/inputs/short_reads/sr.bam --barcodes /homeb/user/AD_data/illumina/illumina-cellranger/HBA8_3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --tes /homeb/user/AD_data/AD-MultiOmics/bulk2cell/results/HBA8_full/tes/raw/tes.tsv --training-cache /homeb/user/AD_data/AD-MultiOmics/bulk2cell/results/HBA8_full/capture_training_prefetch/training.json --tau 0 --out /homeb/user/AD_data/AD-MultiOmics/bulk2cell/results/HBA8_full/quantification_isoseq_prior --workers 8 --max-likelihood-entries 1200000000 "$@"
