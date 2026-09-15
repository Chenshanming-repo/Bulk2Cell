import importlib.util
from pathlib import Path
import json
import pysam

def test_exports_only_primary_called_cell_rna_in_sequenced_orientation(tmp_path):
    script=Path(__file__).parents[2]/'bin/extract_pseudobulk.py'
    assert script.exists(), 'pseudobulk extraction is not implemented'
    spec=importlib.util.spec_from_file_location('pseudobulk',script)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    bam=tmp_path/'reads.bam'
    with pysam.AlignmentFile(bam,'wb',header={'HD':{'VN':'1.6'},'SQ':[{'SN':'chr1','LN':1000}]}) as handle:
        for name,flag,barcode,umi in [('ok',16,'CELL-1','AAA'),('other',0,'OTHER-1','AAA'),('secondary',256,'CELL-1','AAA'),('no_umi',0,'CELL-1','')]:
            read=pysam.AlignedSegment();read.query_name=name;read.query_sequence='AACG';read.query_qualities=pysam.qualitystring_to_array('IIII')
            read.flag=flag;read.reference_id=0;read.reference_start=10;read.cigarstring='4M';read.set_tag('CB',barcode)
            if umi:read.set_tag('UB',umi)
            handle.write(read)
    cells=tmp_path/'barcodes.tsv';cells.write_text('CELL-1\n')
    out=tmp_path/'rna.fq';qc=tmp_path/'qc.json'
    module.extract(bam,cells,out,qc)
    assert out.read_text()=='@ok\nCGTT\n+\nIIII\n'
    assert json.loads(qc.read_text())['accepted_reads']==1
