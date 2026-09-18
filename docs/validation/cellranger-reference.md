# Existing Cell Ranger reference — 2026-09-18

Implemented plan:

1. Add optional `cellranger_reference` with reference-directory structure checks.
2. Share the supplied reference across samples and skip `CELLRANGER_MKREF`; preserve reference building when omitted.
3. Document usage, required genome/annotation inputs and reference compatibility.
4. Verify both workflow routes and invalid inputs with the two-sample smoke test.

Verification with Nextflow 25.10.0:

- `bash tests/workflow/smoke.sh`: passed. Default route completed 24 tasks; supplied-reference route completed 22 tasks with no `CELLRANGER_MKREF` tasks. Both HiFi and FLNC samples reached final stub outputs and staged the same supplied reference, including a directory name containing spaces and compressed reference GTF.
- Missing paths, regular files and incomplete reference directories were rejected before task submission.
- `python -m pytest tests/workflow/test_inputs.py -q`: 16 passed.
- `bash -n tests/workflow/smoke.sh` and `git diff --check`: passed.

These are workflow routing and input checks using artificial fixtures. Real Cell Ranger alignment and biological compatibility of a supplied reference were not tested.
