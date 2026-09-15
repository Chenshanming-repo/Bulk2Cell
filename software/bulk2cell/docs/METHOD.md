# Statistical method and implementation

## Candidate transcripts and molecular evidence

Reference exon GTF and sample Iso-Seq/Pigeon models form a union. Internal coordinates are zero-based half-open. Exact full exon structures are deduplicated within gene, chromosome and strand; reference identifiers take precedence and alias provenance is retained. This deduplication is distinct from 3′ grouping: structurally distinct models remain separate transcript members.

Each molecule is `(cell barcode, corrected UMI, assigned gene)`. All eligible primary sense alignments contribute aligned blocks and exact CIGAR-N junctions. A transcript must contain **all** aligned bases and all splice junctions. Conflicting read evidence can make the entire molecule incompatible. The 3′-most boundary supplies one capture-distance likelihood; PCR/read multiplicity does not multiply independent UMI likelihoods. Deletions contribute no aligned bases and are not splice junctions. The pipeline relies on Cell Ranger's corrected tags and single-gene assignment rather than implementing barcode/UMI correction anew.

## Hybrid abundance

For each gene, normalize Salmon TPM to `w_SR`; absent or all-zero SR mass uses a reported uniform distribution. With LR counts `C` and `N=sum(C)`, use

`w* = (C + tau*w_SR)/(N + tau)`.

At `N=tau=0`, use `w_SR`. `tau` is a nonnegative CLI parameter, shared across genes in this version. Bulk LR counts are support counts potentially biased by amplification/capture; they are not guaranteed independent RNA molecules.

## Capture distribution and TES

Learn a reflected Gaussian KDE on nonnegative **spliced** SR read-to-TES distances from structurally uniquely compatible molecules, pooled across selected genes. Use a seeded reservoir of up to 2,000 molecules and a configurable bandwidth (default 30 nt). When no such molecules exist, record an exponential fallback of scale 300 nt. The model is trained on SR distances to annotation TES; LR endpoint observations are used only for TES marginalization. Calibration uncertainty and bandwidth sensitivity are not estimated automatically.

For transcript k, a weighted endpoint TSV defines normalized `q_k(s)`. Missing data gives a point mass at annotated TES. The molecule likelihood is `L_mk = sum_s q_k(s) f(d(m,s))`, with zero for incompatible structures/negative distance. Computation uses log densities, log-sum-exp, and row normalization to prevent tiny KDE tails from silently dropping compatible molecules. Absolute density scale cancels in each posterior.

## Groups and EM

Start with identical spliced structures within the terminal capture window (default 600 nt), retaining exact TES boundaries. Refine groups when observed likelihood columns differ; upstream evidence or distinct TES distributions can therefore separate members. This is a conservative operational grouping rule, not a proof that different groups are identifiable at every cell's depth. Group IDs are scoped to the catalog, selected evidence and run settings.

Group likelihood is the hybrid-weighted mixture of member likelihoods. Fit each cell/gene using EM with group prior `pi_G = sum_{k in G} w*_k` and regularization `alpha` (`--em-strength`, default 1):

`r_mG = L_mG * p_G / sum_H L_mH * p_H`

`p_G = (sum_m r_mG + alpha*pi_G)/(M + alpha)`.

Initialization is uniform in the interior, permitting data to recover zero-prior components. Stop at max absolute abundance change <1e-7 or the configured iteration limit (default 500); nonconvergence is recorded. Report expected counts `X_cG=sum_m r_mG` without prior pseudomolecules. Split group counts by normalized member hybrid weights. If a group has zero hybrid mass, its within-group split is uniform and its members' zero priors are visible in `isoforms.tsv`.

Group counts themselves may be EM-ambiguous expected counts. Multi-member transcript counts are explicitly labeled `prior_informed_decomposition`; they cannot establish cell-specific changes among members using the shared bulk weights. No posterior credible intervals are claimed. The model does not include read alignment quality scores, internal-priming detection or transcript effective-length correction beyond supplied Salmon TPM.

## Files and responsibilities

| Module | Responsibility |
| --- | --- |
| `models.py` | Immutable transcripts/molecules and transcript structure validation |
| `adapters.py` | GTF/Pigeon/catalog, Salmon, TES, barcodes and BAM readers |
| `tes.py` | Strand-aware empirical endpoints with splice-chain, clipping and terminal-exon support filters |
| `flnc.py` | Exact FLNC selection through PacBio PBI virtual offsets, validated FASTA publication |
| `preparation.py`, `orchestration.py` | Matched annotations/SR inputs, Salmon, validated pbmm2 alignment and source provenance |
| `inference.py` | Compatibility, distances, grouping, KDE, priors and EM |
| `pipeline.py` | Input integration, global distance training, sparse matrices, QC and provenance |
| `evaluation.py` | ID-aligned clustering and abundance metrics, optional UMAP |
| `reporting.py` | File-oriented evaluation with zero-cell exclusion and JSON/TSV reports |
| `cli.py`, `__main__.py` | CLI validation and package entry points |
| `examples/run_hba8.py` | Reproducible bounded real-data launcher |

Unit tests cover known hybrid limits, positive/negative strand splicing, conflicting junctions, TES mixtures, retained identifiability, numerical tails, EM recovery and molecule conservation. Synthetic BAM integration verifies actual CB/UB/GX ingestion, matrix dimensions and overwrite refusal. Evaluation tests check shuffled IDs and split overlap. Real-data support comes from the separately recorded HBA8 benchmark.

## Empirical endpoint evidence

`estimate-tes` keeps raw-read and consensus evidence separate. Raw FLNC observations receive unit weight; consensus support weighting is optional and does not reconstruct within-cluster variation. Eligibility requires primary mapped evidence, MAPQ >=20, at most 20 clipped bases at the three-prime end, the full annotated splice chain within 10 nt, actual aligned-base overlap in the terminal splice exon, and TES within 200 nt of annotation. Insertions, deletions and match/mismatch CIGAR segments do not split an exon; deletion-only overlap provides no aligned-base support. Coordinates remain zero-based interbase boundaries. Missing eligible empirical data retains the annotated point-mass TES.
