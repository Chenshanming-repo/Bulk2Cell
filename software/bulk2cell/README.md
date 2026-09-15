# bulk2cell

<div align="center">
  <img src="../../docs/figures/bulk2cell_workflow.svg" alt="bulk2cell workflow: independent Iso-Seq and Illumina branches converge on single-cell isoform quantification" width="1100">
</div>

A documented Python implementation of **bulk long-read-assisted 10x 3′ isoform quantification**, based on `Bulk_Long-Read-Assisted_Short-Read_Quantification.md`. It combines reference/Iso-Seq structures, optional Salmon abundance and LR TES observations, a smooth short-read capture model, and cell-specific group EM. It is a research implementation; the real-sample benchmark below establishes functionality, not improved biological accuracy.

The existing [Iso-Seq](../isoseq/README.md) and [Cell Ranger](../cellranger/README.md) pipelines remain upstream entry points. The existing [region quantifier](../isoforms_quantification/README.md) measures different, overlapping region targets; its outputs are not substituted for isoform counts.

## Install and test

```bash
python -m venv .venv
.venv/bin/python -m pip install './bulk2cell[test,evaluation]'
.venv/bin/python -m pytest bulk2cell/tests -q
.venv/bin/bulk2cell --help
```

Python 3.10+; NumPy, SciPy and pysam are required. UMAP is optional: install `./bulk2cell[umap]`. A local verified environment is available at `bulk2cell/.venv` on this host. It inherits the existing quantification environment without changing it.

## Salmon installation

Salmon is an external executable, separate from the Python package. The local
installation uses `bulk2cell/.tools/salmon`; reproduce it with:

```bash
conda env create --prefix "$PWD/bulk2cell/.tools/salmon" \
  --file bulk2cell/environment.salmon.yml
bulk2cell/.tools/salmon/bin/salmon --version
```

The installation was verified with Salmon **1.10.3** by building an index and
quantifying 100 synthetic reads, recovering the expected 60/40 transcript counts.
The exact Linux package set is recorded in [salmon-linux-64.lock.txt](salmon-linux-64.lock.txt).
A project-local micromamba installer is available at
`bulk2cell/.tools/micromamba/bin/micromamba` for faster environment creation.
Both commands are available after activating the project Python environment:

```bash
source bulk2cell/.venv/bin/activate
bulk2cell --version
salmon --version
```

Installing Salmon does not generate `quant.sf`. Quantification additionally
requires a matching union-transcriptome index and suitable short-read pseudobulk
input. Supply its output with `bulk2cell quantify --salmon path/to/quant.sf ...`.

## Prepare Iso-Seq benchmark inputs and empirical TES

The reusable preparation commands keep source files read-only and refuse existing outputs:

```bash
bulk2cell prepare --catalog catalog.json.gz --genes GENE1 GENE2 \
  --genome genome.fa --bam possorted_genome_bam.bam \
  --barcodes barcodes.tsv.gz --salmon salmon --out prepared

bulk2cell estimate-tes --isoseq-root sample_isoseq --catalog prepared/catalog.json \
  --out tes_consensus

bulk2cell extract-flnc --isoseq-root sample_isoseq --catalog prepared/catalog.json \
  --out selected.fasta --assignments-json raw_assignments.json
bulk2cell align-flnc --fasta selected.fasta --genome genome.fa \
  --pbmm2 pbmm2 --out flnc.mapped.bam
bulk2cell estimate-tes --mode raw --catalog prepared/catalog.json \
  --bam flnc.mapped.bam --assignments sample_isoseq/05_collapse/sample.collapsed.read_stat.txt \
  --out tes_raw
```

The first TES command discovers the standard `04_align` consensus BAM and collapse group table. Its evidence is labeled `consensus`; optional `--weight-mode support` uses cluster support but cannot reconstruct within-cluster endpoint variation. Raw mode requires an explicitly mapped BAM and read-stat table unless a read-stat file can be discovered; it never treats a clustered consensus BAM as raw evidence. `align-flnc` uses pbmm2 `--preset ISOSEQ --sort --unmapped`, retaining unmapped selected queries for mapping QC. Supply the reference FASTA when an existing `.mmi` was built with other minimizer parameters; the command validates `.mmi` headers and requires ISOSEQ `w=5,k=15`, uncompressed minimizers, and reference sequences. Alignment completion requires a readable indexed BAM and a successful record scan. Each command records input fingerprints, exact external commands where applicable, evidence type, QC, and completion status.

`prepare` always exports union and reference GTF/FASTA, extracts the shared barcode-filtered short-read panel, and runs both targeted Salmon indexes when `--salmon` is supplied. These transcriptome-only indexes support a controlled targeted benchmark; production analysis needs a suitable decoy-aware whole transcriptome.

## Quantify

```bash
bulk2cell quantify \
  --reference reference.gtf.gz \
  --isoseq-gff sample.collapsed.sorted.filtered_lite.gff \
  --classification sample_classification.filtered_lite_classification.txt \
  --bam sample/outs/possorted_genome_bam.bam \
  --barcodes sample/outs/filtered_feature_bc_matrix/barcodes.tsv.gz \
  --salmon salmon/quant.sf --tes lr_tes.tsv \
  --out results/sample
```

The BAM must be coordinate sorted/indexed and carry corrected `CB`, `UB` and a single assigned `GX`. Use `--genes ENSG... ENSG...` or unambiguous gene names for a targeted run. Reference and BAM must use the same genome build/contig names. Missing, ambiguous or antisense gene assignments are excluded and counted; this implementation does not remap reads to novel genes absent from Cell Ranger's reference. To quantify such genes, rerun upstream alignment/assignment against the union reference.

Use `--catalog catalog.json.gz` instead of raw annotations to import an existing regionquant catalog. That adapter inherits the catalog's prior selection/exclusions; raw GTF/Pigeon import can retain novel genes. Exact full structures within the same gene are deduplicated and all original transcript aliases are retained in `run.json`. Distinct upstream structures remain separate members even when grouped. Pigeon import excludes explicit reverse-transcription-switching artifacts, noncanonical models and unsupported structural categories; use upstream Pigeon/SQANTI filtering for internal priming and other quality checks.

`quant.sf` must refer to the **same union annotation**. TPM is normalized within each gene; NumReads is a documented fallback if TPM is absent. Unknown IDs are counted. Aliases from exact-structure deduplication are summed. Without `--salmon`, the SR component is uniform and the fallback is recorded. This is not an inferred Salmon result. To prepare a matching Salmon index, export/use the union `transcripts.gtf`, obtain transcript FASTA with `gffread -g genome.fa -w transcripts.fa transcripts.gtf`, and run Salmon on the short-read pseudobulk using a library configuration appropriate to the data. Salmon is an external dependency. `prepare --salmon EXECUTABLE` runs it explicitly; `quantify` consumes the resulting file.

TES TSV uses genomic **zero-based interbase boundaries**:

```text
transcript_id	tes	count
PB.1.1	123450	8
PB.1.1	123460	2
```

Each row describes observed endpoints of reads assigned to that full isoform; `count` defaults to 1. On +, TES is the alignment's exclusive end; on −, its inclusive start. Supply training LR endpoint observations, not the collapsed annotation endpoint repeated as though it were independent data. Endpoints may move within/extend the terminal exon; a shift into an upstream exon is invalid for that model. Missing TES data uses the annotation boundary and is reported. Use `extract-flnc`, `align-flnc`, and `estimate-tes --mode raw` for endpoint extraction from the original Iso-Seq reads.

## Outputs

Both compressed Matrix Market matrices are **cells × features**, with fractional expected UMI counts. This is the transpose of the standard 10x layout.

| File | Meaning |
| --- | --- |
| `group_counts.mtx.gz`, `groups.tsv` | EM expected group counts and group interpretation |
| `isoform_counts.mtx.gz`, `isoforms.tsv` | Transcript estimates, sources, LR support and hybrid priors |
| `membership.tsv` | Full group membership and fixed within-group decomposition weights |
| `barcodes.tsv` | Matrix row IDs; includes zero-count whitelist cells; header present |
| `gene_qc.tsv` | Assigned/incompatible molecules and EM convergence per gene |
| `transcripts.gtf` | Actual deduplicated union used in this run |
| `distance_training.tsv` | Retained structurally unique SR training distances |
| `run.json` | Completion marker, parameters, aliases, exclusions, input fingerprints, software versions, time and peak RSS |

Existing output directories are refused. An interrupted output without `run.json` is incomplete. Matrix totals are checked against assigned molecules before completion. Large inputs use metadata fingerprints; inputs up to 64 MiB also receive SHA-256. BAM read QC counts are **gene-window fetch events**, which may repeat across overlapping genes; assigned molecules are counted once per `(CB, UB, GX)`.

## Evaluate

```bash
bulk2cell evaluate --run results/sample --out results/sample_evaluation \
  --clusters sample/outs/analysis/clustering/gene_expression_graphclust/clusters.csv \
  --n-clusters 10 --seed 0
```

Library-size normalization, log1p, uncentered sparse TruncatedSVD (PCA-like), deterministic K-means and ARI/NMI/silhouette are provided. Zero-count cells are excluded and reported. Add `--umap` for UMAP coordinates. TSV embeddings and cluster assignments can be plotted in R, Scanpy or Python. Cell Ranger labels are a **comparison reference, not ground truth**. A targeted gene panel is not sufficient for a genome-wide clustering conclusion.

`--original converted_scalpel_directory` compares against independently generated original SCALPEL output. Convert its counts into the documented cell-by-feature `isoform_counts.mtx.gz`, `barcodes.tsv` (`barcode` header), and `isoforms.tsv` (`transcript_id` header) layout first. Cell barcodes are aligned by ID and feature sets may differ. Original SCALPEL is not silently approximated by disabling one model component.

The Python `evaluation.abundance_accuracy()` API compares aligned truth/estimate IDs with MAE, RMSE, Pearson, Spearman, Jensen–Shannon **distance** (square root of divergence) and total variation. For compositional accuracy, supply normalized within-gene proportions. Optional paired `training_read_ids` and `validation_read_ids` reject read overlap. Transcript IDs are expected to overlap. Independent validation requires excluding held-out reads from annotation discovery, counts and TES training, not merely from the final abundance table.

## Reproduce HBA8

```bash
bulk2cell/.venv/bin/python -m pip install --no-deps -e ./bulk2cell
OPENBLAS_NUM_THREADS=2 bulk2cell/.venv/bin/python \
  bulk2cell/examples/run_hba8.py --out bulk2cell/results/HBA8
bulk2cell/.venv/bin/bulk2cell evaluate \
  --run bulk2cell/results/HBA8 --out bulk2cell/results/HBA8_evaluation \
  --clusters /homeb/user/AD_data/illumina/illumina-cellranger/HBA8_3/outs/analysis/clustering/gene_expression_graphclust/clusters.csv \
  --n-clusters 10
```

The example uses 17 genes and all 10,250 filtered barcodes. See the [historical benchmark](docs/BENCHMARK.md) for this earlier uniform-SR example. The [Iso-Seq/Salmon/SCALPEL benchmark](docs/BENCHMARK_ISOSEQ.md) records the completed matched-input comparison and reproduction commands. It reads the existing catalog derived from the provided HBA8 Iso-Seq results. The package can process all catalog genes by omitting `--genes`, but this release retains eligible molecule evidence in memory and dense likelihoods per gene; it is not yet a validated genome-wide scalable replacement for the regionquant sharded runner. `--max-likelihood-entries` bounds individual gene likelihood allocation, not total BAM evidence memory.

## Completed standard SCALPEL comparison

The [two-part HBA8 comparison](results/HBA8_full/standard_comparison/README.md) includes completed standard SCALPEL, previous SCALPEL reference/union and all three Bulk2Cell variants. It reports centered PCA/Leiden agreement with CellRanger gene partitions across two profiles, six resolutions and three seeds, plus fixed cell-type recovery and pooled gene-expression agreement. Rankings reverse with preprocessing; missing standard estimates are explicitly excluded at feature/gene level without imputation.

## Cell-type benchmark

The [HBA8 cell-type benchmark](results/HBA8_full/celltype_benchmark/README.md) compares all three Bulk2Cell variants with the previous SCALPEL reference and union runs using one fixed gene-derived annotation. It includes all 10,250 cells, full clustering sensitivity, per-type recovery, and pooled gene-expression agreement, with independently verified results and downloadable figures. Marker-derived labels measure agreement rather than validated accuracy.

## Method and code map

Read [METHOD.md](docs/METHOD.md) for equations, assumptions and module responsibilities. Every Python module and named function includes a docstring; inline comments explain numerical or coordinate decisions.

Upstream SCALPEL: https://github.com/plasslab/SCALPEL, inspected at `3c31ffc6aa3b7e422623aececfaea32b3aa358ca` (AGPL-3.0). No upstream source files are vendored. This package independently implements the supplied method; its regularization, grouping, KDE and UMI likelihood differ from original SCALPEL. Cite the SCALPEL authors/publication when describing methodological ancestry.

## Empirical TES in HBA8

TES means **transcription end site**, the transcript's 3′ endpoint. Empirical TES
means an endpoint distribution measured from actual long reads assigned to an
isoform: for example, 8 reads ending at position 123450 and 2 at 123460 give
probabilities 0.8 and 0.2. A single annotation coordinate is a point estimate,
not a measured distribution. Here `q_k(s)` models biological/measurement endpoint
variation, while `f(d)` models short-read capture distance to those endpoints.

HBA8 has `05_collapse/HBA8.collapsed.read_stat.txt` mapping original reads to
isoforms, and `02_refine/HBA8.flnc.bam` containing full-length non-concatemer reads.
The inspected `04_align/HBA8.mapped.bam` records use `transcript/...` names and
represent cluster consensus transcripts. Their endpoints can describe consensus
variation, but should not be counted as individual FLNC endpoint observations.
To estimate a read-level distribution, align the original FLNC reads, join their
IDs to the collapse mapping, and collect strand-aware 3′ genomic boundaries for
quality-filtered reads. Truncation, internal priming and uncertain assignments
need filtering; raw endpoints are observations, not automatically verified TES.

## Rename and historical results

The project/package/CLI was renamed from `lrscalpel` to `bulk2cell`. Existing
benchmark files were moved with the project. Their historical `run.json` and
`evaluation.json` records retain the original software name and invocation paths
so the recorded provenance is unchanged; current commands use `bulk2cell`.
