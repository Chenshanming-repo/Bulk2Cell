"""End-to-end CLI checks for standard Iso-Seq output discovery and TES export."""
import json
import subprocess
import sys
from pathlib import Path
import pysam


def test_estimate_tes_cli_uses_standard_isoseq_outputs(tmp_path):
    """A standard consensus BAM/group mapping feeds canonical TES output directly."""
    root=tmp_path/'IsoSeq';(root/'04_align').mkdir(parents=True);(root/'05_collapse').mkdir()
    (root/'05_collapse/s.collapsed.group.txt').write_text('PB.1.1\ttranscript/1\n')
    bam=root/'04_align/s.mapped.bam'
    with pysam.AlignmentFile(bam,'wb',header={'SQ':[{'SN':'chr1','LN':1000}]}) as out:
        r=pysam.AlignedSegment();r.query_name='transcript/1';r.query_sequence='A'*110;r.flag=0
        r.reference_id=0;r.reference_start=100;r.mapping_quality=60;r.cigarstring='110M';out.write(r)
    catalog={'genes':{'g':{'gene_id':'g','gene_name':'G','chrom':'chr1','strand':'+','ref':[],
        'lr':[{'id':'PB.1.1','exons':[[100,200]],'fl':3}]}}}
    cp=tmp_path/'catalog.json';cp.write_text(json.dumps(catalog));dest=tmp_path/'tes'
    cmd=[sys.executable,'-m','bulk2cell','estimate-tes','--isoseq-root',str(root),'--catalog',str(cp),'--out',str(dest)]
    result=subprocess.run(cmd,capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    assert 'PB.1.1\t210\t1.0\tconsensus' in (dest/'tes.tsv').read_text()
    assert json.loads((dest/'run.json').read_text())['evidence']=='consensus'
    assert subprocess.run(cmd,capture_output=True).returncode!=0


def _catalog(path):
    """Write a small catalog containing reference and Pigeon structures."""
    payload={'genes':{'g':{'gene_id':'g','gene_name':'G','chrom':'chr1','strand':'+',
        'ref':[{'id':'R.1','exons':[[100,210]]}],
        'lr':[{'id':'PB.1.1','exons':[[100,200]],'fl':3}]}}}
    path.write_text(json.dumps(payload));return path


def test_estimate_tes_cli_accepts_explicit_raw_evidence(tmp_path):
    """An explicit FLNC alignment produces raw-read evidence with input provenance."""
    from bulk2cell.cli import main
    bam=tmp_path/'raw.bam'
    with pysam.AlignmentFile(bam,'wb',header={'SQ':[{'SN':'chr1','LN':1000}]}) as out:
        r=pysam.AlignedSegment();r.query_name='movie/1/ccs';r.query_sequence='A'*110;r.flag=0
        r.reference_id=0;r.reference_start=100;r.mapping_quality=60;r.cigarstring='110M';out.write(r)
    table=tmp_path/'read_stat.txt';table.write_text('id\tlength\tpbid\nmovie/1/ccs\t110\tPB.1.1\n')
    out=tmp_path/'raw-tes';catalog=_catalog(tmp_path/'catalog.json')
    assert main(['estimate-tes','--mode','raw','--catalog',str(catalog),'--bam',str(bam),
        '--assignments',str(table),'--out',str(out)])==0
    record=json.loads((out/'run.json').read_text())
    assert record['evidence']=='raw_read'
    assert record['inputs']['bam']['path']==str(bam)
    assert record['inputs']['assignments']['path']==str(table)


def test_extract_flnc_cli_uses_pbi_and_publishes_complete_fasta(tmp_path):
    """The extraction command uses an adjacent PBI and verifies every selected query."""
    import gzip,struct
    from bulk2cell.cli import main
    root=tmp_path/'isoseq';(root/'02_refine').mkdir(parents=True);(root/'05_collapse').mkdir()
    name='movie/7/ccs';bam=root/'02_refine/sample.flnc.bam'
    with pysam.AlignmentFile(bam,'wb',header={'HD':{'VN':'1.6'}}) as handle:
        offset=handle.tell();r=pysam.AlignedSegment();r.query_name=name;r.query_sequence='ACGT';r.flag=4;handle.write(r)
    with gzip.open(str(bam)+'.pbi','wb') as pbi:
        pbi.write(struct.pack('<4sIHI18s',b'PBI\x01',0x00040000,0,1,b'\0'*18))
        pbi.write(struct.pack('<iiiifBq',0,0,0,7,1.0,0,offset))
    (root/'05_collapse/sample.collapsed.read_stat.txt').write_text(f'id\tlength\tpbid\n{name}\t4\tPB.1.1\n')
    fasta=tmp_path/'selected.fa';assignments=tmp_path/'assignments.json'
    assert main(['extract-flnc','--isoseq-root',str(root),'--catalog',str(_catalog(tmp_path/'catalog.json')),
        '--out',str(fasta),'--assignments-json',str(assignments)])==0
    assert fasta.read_text()==f'>{name}\nACGT\n'
    record=json.loads(Path(str(fasta)+'.run.json').read_text())
    assert record['mode']=='pbi' and record['selected']==record['found']==1
    assert json.loads(assignments.read_text())=={name:'PB.1.1'}


def test_prepare_cli_builds_real_synthetic_panel_without_salmon(tmp_path):
    """Preparation exports both transcriptomes and matched short reads without Salmon."""
    from bulk2cell.cli import main
    genome=tmp_path/'genome.fa';genome.write_text('>chr1\n'+'A'*500+'\n');pysam.faidx(str(genome))
    bam=tmp_path/'cells.bam'
    with pysam.AlignmentFile(bam,'wb',header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':500}]}) as out:
        r=pysam.AlignedSegment();r.query_name='sr';r.query_sequence='A'*20;r.query_qualities=pysam.qualitystring_to_array('I'*20)
        r.flag=0;r.reference_id=0;r.reference_start=120;r.cigarstring='20M';r.set_tag('CB','cell-1');r.set_tag('UB','umi');r.set_tag('GX','g');out.write(r)
    pysam.index(str(bam));barcodes=tmp_path/'barcodes.tsv';barcodes.write_text('cell-1\n');dest=tmp_path/'prepared'
    assert main(['prepare','--catalog',str(_catalog(tmp_path/'catalog.json')),'--genes','G','--genome',str(genome),
        '--bam',str(bam),'--barcodes',str(barcodes),'--out',str(dest)])==0
    record=json.loads((dest/'run.json').read_text())
    assert record['status']=='complete' and record['salmon']=={} and record['short_read_qc']['accepted_reads']==1
    assert (dest/'union.fa').is_file() and (dest/'reference.gtf').is_file() and (dest/'short_reads/sr.bam.bai').is_file()


def test_align_flnc_cli_rejects_wrong_mmi_and_records_exact_command(tmp_path):
    """Alignment rejects incompatible indexes and records ISOSEQ plus unmapped retention."""
    import struct
    from bulk2cell.cli import main
    fasta=tmp_path/'reads.fa';fasta.write_text('>r\nA\n');bad=tmp_path/'bad.mmi';bad.write_bytes(b'MMI\x02'+struct.pack('<4iI',10,15,14,0,1))
    assert main(['align-flnc','--fasta',str(fasta),'--genome',str(bad),'--out',str(tmp_path/'bad.bam')])==2
    import shlex
    fixture=tmp_path/'fixture.bam'
    with pysam.AlignmentFile(fixture,'wb',header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':100}]}) as handle:
        read=pysam.AlignedSegment();read.query_name='r';read.query_sequence='A';read.flag=0
        read.reference_id=0;read.reference_start=0;read.cigarstring='1M';handle.write(read)
    pysam.index(str(fixture))
    executable=tmp_path/'pbmm2'
    executable.write_text('#!/bin/sh\ncp '+shlex.quote(str(fixture))+' "$4"\ncp '+shlex.quote(str(fixture)+'.bai')+' "$4.bai"\n')
    executable.chmod(0o755)
    genome=tmp_path/'genome.fa';genome.write_text('>chr1\nA\n');out=tmp_path/'mapped.bam'
    assert main(['align-flnc','--fasta',str(fasta),'--genome',str(genome),'--out',str(out),
        '--pbmm2',str(executable),'--threads','3','--sort-threads','1'])==0
    record=json.loads(Path(str(out)+'.run.json').read_text())
    assert record['evidence']=='raw_read'
    assert record['command'][-8:]==['--preset','ISOSEQ','--sort','--unmapped','-j','3','-J','1']
    assert main(['align-flnc','--fasta',str(fasta),'--genome',str(genome),'--out',str(out),'--pbmm2',str(executable)])==2


def test_raw_tes_requires_explicit_bam_and_prepare_rejects_zero_threads(tmp_path):
    """Raw mode never consumes a discovered consensus BAM and invalid threads create nothing."""
    from bulk2cell.cli import main
    root=tmp_path/'isoseq';(root/'04_align').mkdir(parents=True);(root/'05_collapse').mkdir()
    (root/'04_align/sample.mapped.bam').write_bytes(b'not raw')
    (root/'05_collapse/sample.collapsed.read_stat.txt').write_text('id\tlength\tpbid\nr\t1\tPB.1.1\n')
    output=tmp_path/'tes'
    assert main(['estimate-tes','--mode','raw','--isoseq-root',str(root),'--catalog',str(_catalog(tmp_path/'c.json')),'--out',str(output)])==2
    assert not output.exists()
    prepared=tmp_path/'prepared'
    assert main(['prepare','--catalog',str(tmp_path/'c.json'),'--genes','G','--genome','missing',
        '--bam','missing','--barcodes','missing','--threads','0','--out',str(prepared)])==2
    assert not prepared.exists()


def test_align_flnc_rejects_hpc_isoseq_index(tmp_path):
    """An otherwise matching ISOSEQ MMI is rejected when HPC minimizers are enabled."""
    import struct
    from bulk2cell.cli import main
    index=tmp_path/'hpc.mmi';index.write_bytes(b'MMI\x02'+struct.pack('<4iI',5,15,14,1,1))
    fasta=tmp_path/'reads.fa';fasta.write_text('>r\nA\n');out=tmp_path/'mapped.bam'
    assert main(['align-flnc','--fasta',str(fasta),'--genome',str(index),'--out',str(out)])==2
    assert not out.exists()


def test_align_flnc_rejects_empty_successful_tool_output(tmp_path):
    """A successful subprocess cannot publish completion for an invalid BAM."""
    from bulk2cell.cli import main
    fasta=tmp_path/'reads.fa';fasta.write_text('>r\nA\n')
    genome=tmp_path/'genome.fa';genome.write_text('>chr1\nA\n')
    executable=tmp_path/'pbmm2';executable.write_text('#!/bin/sh\ntouch "$4"\n');executable.chmod(0o755)
    out=tmp_path/'invalid.bam'
    assert main(['align-flnc','--fasta',str(fasta),'--genome',str(genome),
        '--out',str(out),'--pbmm2',str(executable)])==2
    assert not Path(str(out)+'.run.json').exists()


def test_align_flnc_preserves_preexisting_index_before_running(tmp_path):
    """A preexisting sidecar is protected even when the BAM itself is absent."""
    from bulk2cell.cli import main
    out=tmp_path/'mapped.bam';index=Path(str(out)+'.bai');index.write_text('preserve')
    executable=tmp_path/'pbmm2';executable.write_text('#!/bin/sh\ntouch "$4"\n');executable.chmod(0o755)
    assert main(['align-flnc','--fasta','missing','--genome','missing.fa','--out',str(out),'--pbmm2',str(executable)])==2
    assert not out.exists()
    assert index.read_text()=='preserve'
