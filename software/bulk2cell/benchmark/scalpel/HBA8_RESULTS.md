# HBA8 SCALPEL benchmark results

Completed 2026-09-12 on the shared analysis server. These runs execute pinned upstream SCALPEL commit `3c31ffc6aa3b7e422623aececfaea32b3aa358ca` through the documented core-stage harness. They are not the complete upstream Nextflow workflow or an independent accuracy validation.

## Runs and exports

| Run | Output | Cells | Representative isoforms | Genes | Isoform nnz | Cell-gene groups | Estimated mass | Runner wall time | Peak child RSS |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| SCALPEL reference, native calibration | `bulk2cell/results/HBA8_scalpel_reference_restart` | 10,250 | 253 | 17 | 22,901 | 12,054 | 15,720.499 | 385.363 s | 860,600 KiB |
| SCALPEL union + reference-calibrated capture distribution | `bulk2cell/results/HBA8_scalpel_union_reference_calibrated_restart` | 10,250 | 555 | 17 | 42,146 | 16,179 | 22,952.084 | 523.230 s | 1,022,312 KiB |

The external `/usr/bin/time -v` measurements were 6:25.49 and 8:43.37, with the same 860,600 and 1,022,312 KiB maximum RSS values. These are shared-server observations rather than isolated performance benchmarks. The runner wall clock starts after Python startup and argument validation; the external time includes them.

Each output contains `isoform_counts.mtx.gz`, `gene_counts.mtx.gz`, `barcodes.tsv`, `isoforms.tsv`, `genes.tsv`, `membership.tsv`, and `matrix_summary.json`. Matrices are cells by features. All 10,250 whitelist cells are retained, including zero rows. `membership.tsv` maps collapsed representatives back to original transcript IDs. A cell-gene group is one barcode/gene combination in upstream support or estimates; it is not a cell count and does not imply convergence.

## Exact execution

Both runs used `LC_ALL=C`, `OPENBLAS_NUM_THREADS=2`, `OMP_NUM_THREADS=2`, and `MKL_NUM_THREADS=2`. The reference command supplied `reference.gtf` and reference Salmon quantification. The union command supplied `union.gtf`, union Salmon quantification, and:

```text
--calibration-probabilities bulk2cell/results/HBA8_scalpel_reference_restart/probabilities.tsv
```

The capture distribution SHA-256 is `a7ba559896fd2c086be5b4dae06174596edfc5a17a4b9b7e48ee7568222de8d6`. It is byte-identical to the distribution from the earlier partial reference run.

The harness ran upstream `gtf_processing.R`, `mapping_filtering.R`, `ip_filtering.R`, `compute_prob.R` for the reference, `fragment_probabilities.R`, and the original parsed `em_algorithm` function. Internal priming was disabled through upstream’s empty-IP branch because no sequence-derived IP table was available. Commands and per-stage elapsed times are in each run’s `commands.json`.

## Calibration and interpretation limits

Only two reference genes had single-supported-transcript calibration molecules: ENSG00000109610 had 23 distinct fragments and ENSG00000264364 had 111. Upstream’s strict `nb.reads < 98th percentile` filter retained only the 23-fragment gene, so native capture calibration is extremely narrow for this targeted panel. The union annotation produced zero qualifying single-supported-transcript calibration genes; its original native-calibration attempt remains preserved in `bulk2cell/results/HBA8_scalpel_union` and correctly produced no distribution. The reported union control therefore uses the same-library reference distribution and carries the exact label “SCALPEL union + reference-calibrated capture distribution.”

Upstream rounds relative abundance to three decimals before the harness multiplies it by positive-probability molecules. Consequently estimated mass need not equal integer supported mass: reference is 15,720.499 versus 15,722 positive molecules (delta -1.501), and union is 22,952.084 versus 22,955 (delta -2.916). Candidate molecule sums were 27,021 and 32,324. The exported gene matrix repeats the same estimated value across the representative rows emitted for a cell/gene, then sparse construction sums them; its total therefore intentionally matches `estimates.tsv` rather than raw molecule totals.

The original EM uses `MAX_IT=30` and abundance-change tolerance 0.01, but does not return iteration counts or convergence status. No per-fit or aggregate convergence claim is made. Salmon TPM+1 within-gene normalization and upstream three-decimal abundance rounding are preserved. Output scaling is the rounded relative abundance multiplied by retained positive-probability per-cell/gene molecule totals, rather than the upstream workflow’s Cell Ranger matrix scaling.

## Validation and publication safety

The focused bridge suite has three tests covering GTF coordinates/stable IDs, zero-cell and collapsed-membership export, and refusal to export incomplete runs. `export.py` now requires `summary.json` before reading run tables or writing any matrix or metadata output. Both real exports were checked for expected dimensions, finite nonnegative values, feature-table lengths, and exact agreement of sparse-matrix totals with `estimates.tsv` and `matrix_summary.json`.
