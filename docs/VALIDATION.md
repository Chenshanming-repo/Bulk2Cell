# Delivery validation — 2026-09-15

Project: `/homeb/user/repos/bulk2cell`.

## Verified

- **124 Python tests passed; 1 skipped**, using this project's Python and Salmon Conda environments. The skip is the optional R/SCALPEL integration, whose separate R installation is not part of this workflow. Fourteen warnings originate from Matplotlib/pyparsing deprecations. See [test log](validation/python-tests.log).
- **24 Nextflow tasks completed** in a two-sample stub run launched from the delivered project. Both the HiFi+primers and FLNC routes reached their own bulk2cell outputs, with correct sample joins. See [smoke log](validation/nextflow-smoke.log) and [trace](validation/nextflow-trace.tsv).
- Real synthetic transcript export and **Salmon decoy-aware indexing/quantification passed**: exact-structure deduplication preserved canonical transcript IDs and 100 reads were assigned. This test uses the same no-length-correction options as the workflow. No sequence-bias correction is enabled, because Salmon rejects it together with `--noLengthCorrection`.
- Input tests cover path resolution, missing files, duplicate sample IDs, unsupported BAM stage/assay, missing HiFi primers and unpaired FASTQ lanes. RNA export tests verify called-cell selection, primary alignment selection, corrected UMIs and restoration of sequencing orientation.
- Nextflow `conda,installed` profile resolves to the correct project-local prefixes.
- All 52 vendored engine/source-support files match the SHA256 snapshot in `software/SOURCE.json`.
- Shell smoke script passes `bash -n`.

## Installed tools

| Tool | Version |
|---|---|
| Nextflow | 25.04.8 |
| Python | 3.11 series |
| Iso-Seq | 4.3.0 |
| Pigeon | 26.2.0 |
| pbmm2 | 1.17.0 |
| Salmon | 1.10.3 |
| Cell Ranger (existing official installation) | 10.1.0 |
| NumPy / SciPy / pysam | 1.26.4 / 1.13.1 / 0.22.1 |

All four Conda prefixes exist under `.conda/`; exact Linux package locks are in `envs/locks/`. Cell Ranger is located at `/homeb/user/cellranger-10.1.0/cellranger` and is configured by `conf/this-host.config`.

## Not yet validated

No full biological run of PacBio processing plus Cell Ranger plus quantification was launched: actual input paths, sample pairing, PacBio processing stage/primers and 10x chemistry still need to be supplied. Stub runs validate scheduling/output contracts, not PacBio processing, Cell Ranger alignment, memory requirements or biological accuracy. Installed executable versions and Cell Ranger options were checked; the complete real-data handoff remains to be validated on a study sample.

The initial workflow uses annotation-derived TES. Individual-read empirical TES estimation is not included. See the README for assay and scientific limitations.
