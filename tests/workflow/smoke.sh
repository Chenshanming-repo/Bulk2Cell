#!/usr/bin/env bash
# Validate channel routing for TWO libraries, including both PacBio entry stages.
set -euo pipefail
project=$(cd "$(dirname "$0")/../.." && pwd)
run=$(mktemp -d "${TMPDIR:-/tmp}/bulk2cell-smoke.XXXXXX")
mkdir -p "$run/fastqs"
touch "$run/hifi.bam" "$run/flnc.bam" "$run/genome.fa" "$run/genes.gtf"
printf '>primer_5p\nACGT\n>primer_3p\nTGCA\n' > "$run/primers.fa"
for sample in A B; do
    for mate in 1 2; do
        touch "$run/fastqs/${sample}_S1_L001_R${mate}_001.fastq.gz"
    done
done
cat > "$run/samples.csv" <<CSV
sample,pacbio_bams,pacbio_stage,primers,fastq_dir,fastq_sample,chemistry
A,hifi.bam,hifi,primers.fa,fastqs,A,SC3Pv3
B,flnc.bam,flnc,,fastqs,B,SC3Pv2
CSV
cd "$run"
nextflow run "$project/main.nf" -profile test -stub-run -ansi-log false \
    --input "$run/samples.csv" --genome "$run/genome.fa" \
    --annotation "$run/genes.gtf" --outdir "$run/results" -work-dir "$run/work"
python - "$run" <<'PY'
import csv,json,sys
from pathlib import Path
root=Path(sys.argv[1])
for sample in ('A','B'):
    result=root/'results'/sample/'bulk2cell/quantification/run.json'
    assert json.loads(result.read_text())['status']=='stub',result
with (root/'results/pipeline_info/trace.txt').open() as handle:
    rows=list(csv.DictReader(handle,delimiter='\t'))
assert len(rows)==24,(len(rows),[r['name'] for r in rows])
assert all(r['status']=='COMPLETED' and r['exit']=='0' for r in rows),rows
# Inspect actual staged mkref inputs, including the upstream producer.
mkrefs=[p.parent for p in (root/'work').glob('*/*/.command.sh') if 'mkdir cr_reference' in p.read_text()]
assert len(mkrefs)==2, mkrefs
for work in mkrefs:
    annotation=work/'annotation.gtf'
    assert annotation.is_file(), f'Cell Ranger must stage the original annotation: {work}'
    producer=annotation.resolve().parent/'.command.sh'
    assert 'reference.fa.fai annotation.gtf' in producer.read_text(), producer
    assert not (work/'union.gtf').exists(), work
print(f'PASS: 24 tasks; both HiFi and FLNC sample paths reached bulk2cell. Artifacts: {root}')
PY

# Reuse one prebuilt reference for both libraries, including paths with spaces.
reference="$run/existing reference"
mkdir -p "$reference/fasta" "$reference/genes" "$reference/star"
touch "$reference/reference.json" "$reference/fasta/genome.fa" "$reference/genes/genes.gtf.gz"
nextflow run "$project/main.nf" -profile test -stub-run -ansi-log false \
    --input "$run/samples.csv" --genome "$run/genome.fa" \
    --annotation "$run/genes.gtf" --cellranger_reference "$reference" \
    --outdir "$run/reused-results" -work-dir "$run/reused-work"
python - "$run" <<'PY'
import csv,json,sys
from pathlib import Path
root=Path(sys.argv[1])
with (root/'reused-results/pipeline_info/trace.txt').open() as handle:
    rows=list(csv.DictReader(handle,delimiter='\t'))
assert len(rows)==22,(len(rows),[r['name'] for r in rows])
assert all(r['status']=='COMPLETED' and r['exit']=='0' for r in rows),rows
assert not any('CELLRANGER_MKREF' in r['name'] for r in rows),rows
assert sum('CELLRANGER_COUNT' in r['name'] for r in rows)==2,rows
for sample in ('A','B'):
    result=root/'reused-results'/sample/'bulk2cell/quantification/run.json'
    assert json.loads(result.read_text())['status']=='stub',result
    assert not (root/'reused-results'/sample/'cellranger_reference').exists()
counts=[p.parent for p in (root/'reused-work').glob('*/*/.command.sh')
        if 'mkdir -p count/outs/filtered_feature_bc_matrix' in p.read_text()]
assert len(counts)==2,counts
for work in counts:
    assert (work/'existing reference').resolve()==root/'existing reference',work
print(f'PASS: 22 tasks; both samples reuse the supplied reference without mkref. Artifacts: {root}')
PY

# Invalid inputs fail before any processes are submitted.
mkdir "$run/incomplete-reference"
for invalid in "$run/missing-reference" "$run/genome.fa" "$run/incomplete-reference"; do
    if nextflow run "$project/main.nf" -profile test -stub-run -ansi-log false \
        --input "$run/samples.csv" --genome "$run/genome.fa" \
        --annotation "$run/genes.gtf" --cellranger_reference "$invalid" \
        --outdir "$run/invalid-results" -work-dir "$run/invalid-work" \
        > "$run/invalid.log" 2>&1; then
        echo "FAIL: accepted invalid reference $invalid" >&2
        exit 1
    fi
    if ! grep -Eq 'must be a Cell Ranger reference directory|does not exist|No such file or directory' "$run/invalid.log"; then
        cat "$run/invalid.log" >&2
        exit 1
    fi
    if grep -q 'Submitted process' "$run/invalid.log"; then
        cat "$run/invalid.log" >&2
        exit 1
    fi
done
echo 'PASS: missing, non-directory and incomplete references rejected before task submission.'
