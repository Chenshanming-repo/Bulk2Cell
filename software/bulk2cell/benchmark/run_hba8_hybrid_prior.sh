#!/usr/bin/env bash
# Active HBA8 hybrid prior: LR plus Salmon, tau10;20effective fit workers.
set -euo pipefail
cd /homeb/user/AD_data/AD-MultiOmics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
exec taskset --cpu-list 0-29 bulk2cell/.venv/bin/python bulk2cell/benchmark/run_full_fit_workers.py --fit-workers 20 -- quantify-full --catalog bulk2cell/results/HBA8_full/inputs/catalog.json.gz --bam bulk2cell/results/HBA8_full/inputs/short_reads/sr.bam --barcodes /homeb/user/AD_data/illumina/illumina-cellranger/HBA8_3/outs/filtered_feature_bc_matrix/barcodes.tsv.gz --salmon bulk2cell/results/HBA8_full/inputs/union_salmon/quant/quant.sf --tes bulk2cell/results/HBA8_full/tes/raw/tes.tsv --training-cache bulk2cell/results/HBA8_full/capture_training_prefetch/training.json --out bulk2cell/results/HBA8_full/quantification --workers 8 --max-likelihood-entries 1200000000 --tau 10 "$@"
