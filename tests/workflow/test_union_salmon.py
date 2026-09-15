"""Small real reference-export/decoy-index/Salmon test; no downloaded data."""
import importlib.util
import json
import os
from pathlib import Path
import random
import shutil
import subprocess
import csv
import pysam
import pytest


def test_union_ids_survive_real_decoy_aware_salmon(tmp_path):
    salmon=os.environ.get('SALMON') or shutil.which('salmon')
    if not salmon: pytest.skip('Set SALMON or put salmon on PATH for real Salmon integration')
    root=Path(__file__).parents[2]
    spec=importlib.util.spec_from_file_location('support',root/'bin/workflow_support.py')
    support=importlib.util.module_from_spec(spec);spec.loader.exec_module(support)
    rng=random.Random(41); sequence=''.join(rng.choices('ACGT',k=3000))
    genome=tmp_path/'genome.fa';genome.write_text('>chr1\n'+sequence+'\n');pysam.faidx(str(genome))
    ref=tmp_path/'ref.gtf'
    ref.write_text('chr1\ttest\texon\t1\t800\t.\t+\t.\tgene_id "G1"; transcript_id "T1";\n'
                   'chr1\ttest\texon\t1001\t1800\t.\t+\t.\tgene_id "G2"; transcript_id "T2";\n')
    gff=tmp_path/'lr.gff';gff.write_text('chr1\ttest\texon\t1\t800\t.\t+\t.\tgene_id "PB.1"; transcript_id "PB.1.1";\n')
    classification=tmp_path/'classification.txt'
    classification.write_text('isoform\tstructural_category\tassociated_gene\tfl_assoc\nPB.1.1\tfull-splice_match\tG1\t10\n')
    support.export_union(ref,gff,classification,genome,tmp_path/'union')
    manifest=json.loads((tmp_path/'union.json').read_text())
    assert manifest['annotation_qc']['aliases']['PB.1.1']=='T1'
    fasta=(tmp_path/'union.fa').read_text()
    assert '>T1\n' in fasta and '>T2\n' in fasta and '>PB.1.1' not in fasta
    (tmp_path/'gentrome.fa').write_text(fasta+genome.read_text())
    (tmp_path/'decoys.txt').write_text('chr1\n')
    fastq=tmp_path/'rna.fastq'
    with fastq.open('w') as handle:
        for i in range(100):
            start=(0 if i<60 else 1000)+rng.randrange(50,650)
            handle.write(f'@read{i}\n{sequence[start:start+75]}\n+\n'+75*'I'+'\n')
    with (tmp_path/'salmon.log').open('w') as log:
        subprocess.run([salmon,'index','-t',str(tmp_path/'gentrome.fa'),'-d',str(tmp_path/'decoys.txt'),'-i',str(tmp_path/'index'),'--keepDuplicates','-p','2'],stdout=log,stderr=log,check=True)
        subprocess.run([salmon,'quant','-i',str(tmp_path/'index'),'-l','U','-r',str(fastq),'--validateMappings','--noLengthCorrection','--fldMean','200','--fldSD','80','-p','2','-o',str(tmp_path/'quant')],stdout=log,stderr=log,check=True)
    with (tmp_path/'quant/quant.sf').open() as handle: rows=list(csv.DictReader(handle,delimiter='\t'))
    assert {row['Name'] for row in rows}=={'T1','T2'}
    assert sum(float(row['NumReads']) for row in rows)==pytest.approx(100,abs=1)
