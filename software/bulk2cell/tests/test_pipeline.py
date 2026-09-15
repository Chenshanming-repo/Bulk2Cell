"""End-to-end synthetic CLI tests using a real indexed tagged BAM."""
import json
import subprocess
import sys
import gzip
import pysam
from scipy.io import mmread


def test_cli_conserves_umis_and_refuses_overwrite(tmp_path):
    """An indexed BAM and catalog produce documented sparse cell-by-feature counts."""
    catalog={'genes':{'g':{'gene_id':'g','gene_name':'G','chrom':'chr1','strand':'+',
        'ref':[{'id':'a','exons':[[100,400]]}],
        'lr':[{'id':'b','exons':[[200,400]],'fl':3,'strict':True}]}}}
    cp=tmp_path/'catalog.json'; cp.write_text(json.dumps(catalog))
    bp=tmp_path/'cells.tsv'; bp.write_text('cell-1\nempty-1\n')
    bam=tmp_path/'reads.bam'
    with pysam.AlignmentFile(bam,'wb',header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000}]}) as out:
        for i in range(3):
            r=pysam.AlignedSegment(); r.query_name=f'r{i}'; r.query_sequence='A'*20
            r.flag=0; r.reference_id=0; r.reference_start=350; r.mapping_quality=255
            r.cigarstring='20M'; r.set_tag('CB','cell-1'); r.set_tag('UB',f'u{i}'); r.set_tag('GX','g'); out.write(r)
    pysam.index(str(bam))
    dest=tmp_path/'output'
    args=[sys.executable,'-m','bulk2cell','quantify','--catalog',str(cp),'--bam',str(bam),'--barcodes',str(bp),'--out',str(dest),'--window','50','--tau','0']
    run=subprocess.run(args,capture_output=True,text=True)
    assert run.returncode==0,run.stderr
    assert mmread(dest/'isoform_counts.mtx.gz').sum()==3
    assert mmread(dest/'group_counts.mtx.gz').shape==(2,1)
    manifest=json.loads((dest/'run.json').read_text())
    assert manifest['qc']['assigned_molecules']==3
    assert manifest['parameters']['salmon'] is None
    assert subprocess.run(args,capture_output=True).returncode!=0
