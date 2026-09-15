#!/usr/bin/env python3
"""Validate paired libraries and export the engine's exact transcript union."""
import argparse
import csv
import json
import re
from pathlib import Path


def validate_samples(samplesheet):
    """Resolve local paths relative to CSV; refuse ambiguous sample/library inputs."""
    sheet = Path(samplesheet).resolve()
    required = {'sample','pacbio_bams','pacbio_stage','primers','fastq_dir','fastq_sample','chemistry'}
    rows, seen = [], set()
    with sheet.open(newline='') as handle:
        reader = csv.DictReader(handle)
        if not required.issubset(reader.fieldnames or []):
            raise ValueError('samplesheet requires columns: '+','.join(sorted(required)))
        for number, raw in enumerate(reader, 2):
            row = {key: (raw.get(key) or '').strip() for key in required}
            sample = row['sample']
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}', sample):
                raise ValueError(f'row {number}: invalid sample identifier')
            if sample in seen: raise ValueError(f'duplicate sample: {sample}')
            seen.add(sample)
            if row['pacbio_stage'] not in {'hifi','flnc'}:
                raise ValueError(f'{sample}: stage must be hifi or flnc; convert subreads to HiFi first')
            if row['chemistry'] not in {'SC3Pv2','SC3Pv3','SC3Pv4'}:
                raise ValueError(f'{sample}: chemistry must be SC3Pv2, SC3Pv3 or SC3Pv4 (3-prime GEX)')
            if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]*', row['fastq_sample']):
                raise ValueError(f'{sample}: invalid fastq_sample')
            def resolve(value, directory=False):
                path = Path(value).expanduser()
                if not path.is_absolute(): path = sheet.parent/path
                path = path.resolve()
                if not (path.is_dir() if directory else path.is_file()):
                    raise ValueError(f'{sample}: input not found: {path}')
                # Paths enter shell commands quoted; reject embedded shell syntax/newlines.
                if any(c in str(path) for c in "'\"`$\n\r\t"):
                    raise ValueError(f'{sample}: unsupported characters in path: {path}')
                return str(path)
            bams = row['pacbio_bams'].split(';')
            if any(not bam.strip() or not bam.strip().endswith('.bam') for bam in bams):
                raise ValueError(f'{sample}: pacbio_bams must contain semicolon-separated BAM paths')
            row['pacbio_bams'] = [resolve(bam.strip()) for bam in bams]
            if len(set(row['pacbio_bams'])) != len(bams): raise ValueError(f'{sample}: duplicate BAM')
            if row['pacbio_stage']=='hifi' and not row['primers']:
                raise ValueError(f'{sample}: hifi requires primers FASTA (one 5p/3p pair)')
            if row['primers']:
                row['primers'] = resolve(row['primers'])
                headers = [line[1:].strip() for line in Path(row['primers']).read_text().splitlines() if line.startswith('>')]
                if len(headers)!=2 or not any(h.endswith('5p') for h in headers) or not any(h.endswith('3p') for h in headers):
                    raise ValueError(f'{sample}: primers must contain one header ending 5p and one ending 3p')
            row['fastq_dir'] = resolve(row['fastq_dir'], directory=True)
            pattern = re.compile(re.escape(row['fastq_sample'])+r'_S\d+(?:_L\d{3})?_R([12])_\d{3}\.fastq\.gz$')
            reads = {1:set(),2:set()}
            for path in Path(row['fastq_dir']).glob('*.fastq.gz'):
                match = pattern.fullmatch(path.name)
                if match: reads[int(match[1])].add(re.sub(r'_R[12]_', '_RX_', path.name))
            if not reads[1] and not reads[2]: raise ValueError(f'{sample}: no matching 10x FASTQ files')
            if not reads[1] or reads[1]!=reads[2]: raise ValueError(f'{sample}: FASTQ R1/R2 lanes must be paired')
            rows.append(row)
    if not rows: raise ValueError('samplesheet is empty')
    return rows


def export_union(reference, gff, classification, genome, prefix):
    """Use the same import/deduplication code as final quantification."""
    from bulk2cell.adapters import load_transcript_catalog
    from bulk2cell.preparation import export_transcriptome
    transcripts, qc = load_transcript_catalog(reference, gff, classification)
    result = export_transcriptome(transcripts, genome, prefix)
    Path(str(prefix)+'.json').write_text(json.dumps({'export':result,'annotation_qc':qc}, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command',required=True)
    validate = commands.add_parser('validate')
    validate.add_argument('samplesheet'); validate.add_argument('output')
    union = commands.add_parser('union')
    for name in ('reference','gff','classification','genome','prefix'): union.add_argument('--'+name,required=True)
    args = parser.parse_args()
    if args.command=='validate':
        Path(args.output).write_text(json.dumps(validate_samples(args.samplesheet),indent=2))
    else: export_union(args.reference,args.gff,args.classification,args.genome,args.prefix)

if __name__=='__main__': main()
