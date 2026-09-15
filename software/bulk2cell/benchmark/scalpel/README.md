# Pinned upstream SCALPEL benchmark

This harness executes [SCALPEL](https://github.com/plasslab/SCALPEL) source at
`3c31ffc6aa3b7e422623aececfaea32b3aa358ca` (the actual clone origin is recorded by
`git -C bulk2cell/.tools/scalpel/upstream remote -v`). Upstream is AGPL-3.0;
its source and LICENSE remain in `.tools/scalpel/upstream`, outside the Python
package. This is a scoped upstream-algorithm benchmark, not the complete Nextflow
workflow or an independent accuracy validation.

## Environment and execution

Install the core R dependencies into an isolated prefix:

```bash
MAMBA_ROOT_PREFIX=$PWD/bulk2cell/.tools/mamba-root \
  bulk2cell/.tools/micromamba/bin/micromamba --no-rc create -y \
  -p "$PWD/bulk2cell/.tools/scalpel/env" --override-channels \
  -c https://conda.anaconda.org/conda-forge \
  -c https://conda.anaconda.org/bioconda \
  r-base=4.3 r-argparse r-data.table r-dplyr r-tidyr r-stringr \
  r-ggplot2 bioconductor-genomicranges r-scales

bulk2cell/.venv/bin/python bulk2cell/benchmark/scalpel/run.py \
  --upstream bulk2cell/.tools/scalpel/upstream \
  --rscript bulk2cell/.tools/scalpel/env/bin/Rscript \
  --gtf bulk2cell/results/HBA8_benchmark_inputs/reference.gtf \
  --quant bulk2cell/results/HBA8_benchmark_inputs/reference_salmon/quant.sf \
  --bam bulk2cell/results/HBA8_benchmark_inputs/sr.bam \
  --barcodes bulk2cell/results/HBA8_benchmark_inputs/barcodes.tsv \
  --out bulk2cell/results/HBA8_scalpel_reference
```

Repeat with union GTF/Salmon and a fresh `HBA8_scalpel_union` output for the
matched-annotation control. Existing output directories are refused.



## Full-sample bounded execution

For `bulk2cell/results/HBA8_full/inputs`, use `reference.gtf`,
`short_reads/sr.bam`, `short_reads/barcodes.tsv`, and
`reference_salmon/quant/quant.sf`; use the corresponding union GTF and
`union_salmon/quant/quant.sf` for the union run. Annotation and BED generation
stream into the original chromosome boundaries through a bounded handle pool.
Every selected chromosome alignment is retained, including records that do not
overlap an exon. This preserves upstream `CB::UB` family rejection,
cross-gene UMI collisions, chromosome-range filtering and spliced/unspliced
consistency. The pinned preprocessing stages run one chromosome at a time.
The original EM runs separately for each chromosome after fragment probability
calculation; the runner rejects any gene spanning chromosomes before using this
independence. `compute_prob.R` still runs exactly once on the deterministic
concatenation of all chromosome calibration records, keeping native calibration
global across the full sample.

For chromosomes whose overlap expansion is unsafe, opt into exact whole-cell
batches:

```bash
# Start serially; choose the final size only after measuring a MALAT1 batch.
... run.py ... --boundary cell-batch --cells-per-batch 512 \
  --prepared-reads bulk2cell/results/HBA8_full/shared_reads --r-workers 1
```

A completed shared-read bundle is accepted only after its atomic certificate,
BAM/barcode/consumer fingerprints, selection rules, cross-CB read-ID result and
every required contig BED fingerprint validate. This reuses the same immutable
BED source for reference and union, avoiding another BAM scan and GNU sort. Use
`--stop-after-prepare` for the required measured mapping probe, then remove that
flag and add `--resume` to continue the same source-bound run directory. A
`PREPARED_ONLY.json` marker explicitly states that scientific output is incomplete.

Every CB is assigned to one batch and each batch uses the complete chromosome
exon table. A disk-backed GNU sort rejects a read ID occurring under multiple
CBs, using the text before the first `/` exactly as upstream does. Batch-local
`unique.reads` are discarded: observed gene/transcript pairs are reconciled
globally, corrected local calibration rows are emitted, and the unmodified
`compute_prob.R` runs once. `--r-workers` defaults to one. Increase it only after
measuring a representative high-complexity batch and checking host memory; bounded
workers have exact barriers after chromosome annotation, mapping/IP, global
concordance and calibration, fragment weighting, and EM. Worker failure cancels
queued work and terminates owned R process groups. `mapped.rds` is transient; it
is removed only after filtered RDS and read-ID outputs have a completion
manifest. The fixed upstream distance diagnostic is linked to `/dev/null` and
that disposal is recorded because it is not consumed downstream.

The runner also preserves an upstream concordance bug rather than silently
fixing it. Upstream evaluates `ftrs %in% trs.todel` where `trs.todel` is a
one-column data table: it deletes a candidate only when the chromosome-global
table has exactly one row; tables with zero or multiple rows match nothing. An
instrumented, source-derived first pass records local candidates. If their
chromosome-global distinct count exceeds one, only batches with a local
singleton are rerun through a source-derived no-op bridge. The generated scripts,
upstream source, harness and provenance are bound to resume identity. The old17
gate reproduced native probabilities byte-for-byte and reproduced all 18,552
support and 52,080 estimate data rows exactly after canonical sorting.

Run native calibration independently for reference and union. If native union
calibration fails because upstream finds no eligible single-model calibration
genes, preserve that failure and start a separately labelled
`--calibration-probabilities REFERENCE_RUN/probabilities.tsv` control. Never
replace a successful native union calibration with that control.

Use `--resume` for an interrupted output directory. Completed stages have
`*.complete.json` manifests bound to SHA-256 input contents, parameters,
harness code and the pinned upstream source manifest. Resume takes an exclusive
writer lock and rejects changed identities, missing outputs, or already exported
matrices. It revokes stale completion markers before merging final tables.
Every manifest fingerprints its declared outputs, so present but changed artifacts
are rejected. `invocations.json` accumulates completed invocation wall time; an
unclosed interrupted attempt is reported as duration-unknown and makes
`wall_seconds_complete` false. `summary.json` separates actual bounded-stage spans
from summed command time and reports adopted preparation time explicitly.
Empty-evidence chromosomes are explicit and retain their exon models.

A reviewed runner upgrade never silently rebinds an existing prepared directory.
`adopt_prepared.py` requires the old runner recorded by the legacy manifest,
byte-identical preparation functions, exact old GTF/quant/BAM/barcodes and batch
size, the original certificate lineage, and hashes of every prepared output. It
archives the old manifest outside the active `*.complete.json` glob, writes the
new fingerprinted manifest, and publishes `PREPARE_ADOPTION_COMPLETE.json` last.
The runner refuses an incomplete or stale adoption transaction.

## What is original and what is adapted

* Unmodified upstream `gtf_processing.R`, `mapping_filtering.R`,
  `ip_filtering.R`, `compute_prob.R`, and `fragment_probabilities.R` are run with
  their ordinary CLI. Every command, working directory, exit status, and runtime
  appears in `commands.json`; all stage outputs remain available for inspection.
* `em_batch.R` parses and evaluates only the original `MAX_IT` assignment and
  `em_algorithm` definition from the pinned source, then invokes that same function
  for each cell and gene. No EM equations are reimplemented or source text changed.
  As upstream, relative abundance is rounded to three decimals. The 30-iteration
  cap and 0.01 abundance-change tolerance are unchanged. Upstream does not expose
  iteration counts; convergence per fit is therefore not claimed.
* Truncation distance is 600 nt, isoform end-collapse distance 30 nt, probability
  bin width 20 nt, training-gene count percentile `98%`, and IP distance 60 nt,
  all upstream defaults. No artificial distribution substitutes for `compute_prob`.
* Python emits the nine-column exon table otherwise made by rtracklayer, preserving
  one-based closed GTF coordinates. It applies upstream `(Salmon TPM + 1)` divided
  by the gene total for the single supplied Salmon sample. Safe reversible aliases
  avoid underscore splitting in SCALPEL collapsed-isoform IDs; aliases.json maps
  them to the original gene/transcript IDs. No biological isoform identity changes.
* Input is a shared selected-gene, whitelist-cell, primary, single-GX BAM prepared
  for both methods, not the complete Cell Ranger BAM. Python extracts CB/UB tags
  by name and splits CIGAR at N, in place of positional-tag `bam2bed` shell code.
  No additional samtools markdup runs; original mapping_filtering positional
  fragment deduplication still runs. BAM CB/UB errors and excluded records are
  counted in summary.json. This preprocessing restricts the comparison's scope.
* A header-only internal-priming table selects upstream's existing empty-IP branch.
  Internal-priming sequence filtering is disabled in both annotation comparisons;
  the benchmark does not evaluate that SCALPEL capability.
* Upstream's gene-expression scaling from a Cell Ranger matrix and Seurat object
  construction are replaced by transparent output tables: `estimated_count` is
  the original rounded isoform abundance times the number of distinct cell/gene
  UMIs with positive model probability. `support.tsv` records candidate and
  positive molecules. Counts may have rounding drift. Zero-probability groups
  discarded by upstream have no estimates. Cell Ranger gene calls are not truth.
* `estimates.tsv` carries stable transcript/gene IDs and barcodes. It describes
  collapsed representative isoforms, not automatically all original catalog
  transcripts. Compare gene totals and callability before comparing isoforms.

`summary.json` records wall time, bounded parallel spans, summed command time,
maximum child-process RSS (Linux KiB; largest single child, not a sum),
selected-record QC, parameters/limitations and inputs.
The benchmark reuses LR annotations and SR-derived Salmon pseudobulk, so these are
shared inputs, not independent held-out evidence.

Export completed outputs to sparse matrices only after `summary.json` exists; interrupted runs are refused before any export file is written. All 10,250 whitelist cells are retained,
including zero rows):

```bash
bulk2cell/.venv/bin/python bulk2cell/benchmark/scalpel/export.py \
  bulk2cell/results/HBA8_scalpel_reference \
  bulk2cell/results/HBA8_benchmark_inputs/barcodes.tsv
```

`isoform_counts.mtx.gz` and `gene_counts.mtx.gz` have cells as rows, in headered
`barcodes.tsv` order, and features as columns, in `isoforms.tsv` / `genes.tsv`
order. `membership.tsv` maps each retained collapsed representative to original
transcript IDs; representatives may share members when upstream Salmon ties retain
multiple representatives. `matrix_summary.json` records dimensions and totals.
The environment's exact package URLs are in `environment-linux-64.explicit.txt`.

Pinned source hashes in `upstream-sha256.json` are checked before execution.
Validation completed: three Python bridge tests (GTF coordinates/stable IDs, sparse
zero-cell alignment/collapsed memberships, incomplete-run export refusal); an actual six-stage upstream synthetic
smoke; and an ambiguous two-isoform example whose original EM CLI and shared-process
EM both returned exactly 0.754 / 0.246 after upstream rounding.
