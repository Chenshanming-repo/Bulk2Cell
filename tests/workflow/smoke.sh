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
