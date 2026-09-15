#!/usr/bin/env python3
"""Stream called-cell RNA reads into Salmon FASTQ, preserving sequencing orientation."""
import argparse
from collections import Counter
import json
from pathlib import Path
import pysam
from bulk2cell.adapters import load_barcodes


def extract(bam_path, barcode_path, fastq_path, qc_path):
    """Keep one primary record per RNA read; PCR copies remain abundance evidence."""
    cells=load_barcodes(barcode_path)
    if not cells: raise ValueError('Cell Ranger returned no called cells')
    qc=Counter()
    with pysam.AlignmentFile(str(bam_path),'rb') as bam, open(fastq_path,'w') as output:
        for read in bam.fetch(until_eof=True):
            qc['input_records']+=1
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
                qc['nonprimary_or_unmapped']+=1;continue
            if not read.has_tag('CB') or read.get_tag('CB') not in cells:
                qc['barcode_excluded']+=1;continue
            if not read.has_tag('UB') or not read.get_tag('UB'):
                qc['missing_umi']+=1;continue
            # Standard 10x 3-prime BAM contains the RNA mate only. Refuse a paired
            # representation to avoid exporting barcode sequence or both mates.
            if read.is_paired: raise ValueError('Expected single RNA-mate records in a standard 10x 3-prime Cell Ranger BAM')
            sequence=read.get_forward_sequence(); qualities=read.get_forward_qualities()
            if not sequence or qualities is None:
                qc['missing_sequence_or_quality']+=1;continue
            output.write(f'@{read.query_name}\n{sequence}\n+\n{pysam.array_to_qualitystring(qualities)}\n')
            qc['accepted_reads']+=1
    Path(qc_path).write_text(json.dumps(dict(qc),indent=2))
    if not qc['accepted_reads']: raise ValueError('No called-cell RNA reads available for Salmon')

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ('bam','barcodes','fastq','qc'): parser.add_argument(name)
    args=parser.parse_args();extract(args.bam,args.barcodes,args.fastq,args.qc)
