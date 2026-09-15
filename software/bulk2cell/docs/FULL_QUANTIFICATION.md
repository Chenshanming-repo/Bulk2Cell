# Bounded full-sample quantification

`quantify-full` runs the empirical-TES method, annotation-TES ablation, and exact SR-only ablation together. It retains every selected catalog transcript and whitelist cell, including zero-count features. It writes the standard matrix/TSV interfaces under `empirical_tes/`, `annotation_tes/`, and `sr_only/`.

```bash
OPENBLAS_NUM_THREADS=1 bulk2cell quantify-full \
  --catalog inputs/catalog.json.gz --bam inputs/short_reads/sr.bam \
  --barcodes inputs/short_reads/barcodes.tsv \
  --salmon inputs/union_salmon/quant/quant.sf --tes empirical/tes.tsv \
  --out quantification --workers 8 --max-likelihood-entries 100000000
```

The empirical variant falls back to annotation endpoints for transcripts without TES observations, as in `quantify`. Without `--tes`, empirical and annotation variants are equivalent. The SR-only variant uses the same union structures, Salmon abundance, annotation endpoints, and capture model, with LR abundance set to zero.

The first pass fetches one gene at a time from the indexed BAM. It applies the original structural-uniqueness rule and the same 2,000-distance reservoir, in catalog gene order and sorted molecule order. Training is global across all selected genes; it is never fitted separately for shards. The second pass refetches each gene and fits it with that shared model. Repeated evidence shapes reuse the original scalar likelihood kernel; equal-row cell EM batches preserve each cell's arithmetic and independent stopping iteration. Annotation and SR-only fits share the annotation likelihood.

At most eight worker processes run simultaneously, with one BLAS thread each and at most 128 outstanding training tasks. Training results are consumed in their original order, preserving reservoir sampling exactly. Fitting retains at most twice the worker count of outstanding tasks because its results are larger. `--max-likelihood-entries` bounds each gene's molecule-by-transcript matrix; exceeding it raises an error naming the gene and dimensions during training, before the fit pass. Several arrays and read evidence coexist, so this parameter is an allocation guard, not a complete RSS limit. For 100 million entries, one float64 matrix is 800 MB; budget multiple such arrays per worker. Sparse per-gene results are checkpointed in SQLite and merged by sorted gene IDs without whole-sample dense matrices or coordinate lists.

To overlap capture training with upstream Salmon/TES work:

```bash
OPENBLAS_NUM_THREADS=1 bulk2cell quantify-full \
  --catalog inputs/catalog.json.gz --bam inputs/short_reads/sr.bam \
  --barcodes inputs/short_reads/barcodes.tsv \
  --out capture_training --train-only --workers 8 \
  --max-likelihood-entries 100000000

# Later, add to the ordinary quantify-full command:
# --training-cache capture_training/training.json
```

The training cache is independently bound to evidence fingerprints, selected genes/order, seed, and source code. It does not depend on Salmon, TES, LR/SR prior strength, or EM settings. The later fit run keeps its complete identity, including priors, TES, the training-cache fingerprint, and inference parameters. Large inputs use resolved path, size, and nanosecond mtime; inputs up to 64 MiB also receive SHA-256 hashes. Large-source fingerprints are metadata checks, not cryptographic content verification. Code changes invalidate reuse.

Pass `--resume` with the same output directory after interruption. The input/code/scientific identity must match. Worker count and the likelihood allocation budget may change because they are execution controls. A lock prevents concurrent writers. SQLite uses WAL with NORMAL synchronization to avoid a disk flush per gene. Transactions remain consistent; abrupt power loss may discard recent checkpoints, which are then recomputed. Completed gene checkpoints are reused; missing fits are recomputed. Completion manifests are removed before resumed work and published after output merging and molecule-conservation checks. No failed gene is converted to zeros. Keep checkpoint files local and trusted: their compressed Python payloads are not a format for importing untrusted data.

Runtime is explicitly scoped to the shared three-variant runner, not three independent method timings. Manifests include cumulative observed checkpoint wall time, summed training-worker time, shared extraction/annotation-likelihood time, per-variant worker fit time, and separately reported parent and worker lifetime peak RSS. `elapsed_seconds` includes the current runner's merge; externally cached training elapsed time is reported separately. Time after the last checkpoint before an abrupt interruption can be unrecorded. Worker lifetime peaks are not simultaneous aggregate process-tree RSS. These fields must not be presented as three independent wall-time/RSS measurements against SCALPEL.

## Completion-order fitting launcher

For genes with very different runtimes, the original ordered fit queue can leave workers idle behind one slow gene. The repository launcher `bulk2cell/benchmark/run_full_completion_order.py` consumes completed fits immediately, with the same eight-worker maximum and at most sixteen outstanding tasks. It forwards training to the original ordered scheduler. Per-gene scientific calculations are unchanged, and final outputs still merge in sorted gene order.

From the repository root, use `OPENBLAS_NUM_THREADS=1 bulk2cell/.venv/bin/python bulk2cell/benchmark/run_full_completion_order.py -- quantify-full` followed by the same fitting arguments shown above. For later resumes, use this same launcher and add `--resume`. The adjacent `OUT.scheduler.json` binds the launcher, frozen modules, scientific identity, worker policy and any transition record. A changed launcher, module, policy or transition record is rejected; the original runner's identity checks and writer lock remain active.

Adopting an existing ordinary-run checkpoint additionally requires `--adopt-existing --transition-record FILE` before `--`. Stop its original writer first, create a consistent SQLite backup with the SQLite backup API, and record its verification, original execution, transition time, and uncommitted-work accounting in that JSON record. Preserve the record and backup. The launcher does not create or verify that backup for the caller. Its interrupt handler terminates only its own workers; committed genes survive and uncommitted genes recompute.

The HBA8 transition is recorded in `results/HBA8_full/scheduler_transition.json`. All 524 original committed payloads were independently verified unchanged after resume. The real 17-gene equivalence check at `results/HBA8_full_runner_completion_validation/scheduler_equivalence.json` confirms exact equality of all scientific checkpoint payloads, 21 annotation/QC tables, and six sparse matrices. The observed pilot timing improvement is not a controlled speedup measurement because system load differed. Full-run timing must disclose the original and resumed execution attempts.
