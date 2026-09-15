"""Tests for matched benchmark inputs, transcript coordinates and Salmon preparation."""
from pathlib import Path
import gzip
import json
import pysam
import pytest
from bulk2cell.models import Transcript
from bulk2cell.preparation import export_transcriptome, subset_catalog, extract_short_reads, extract_short_reads_streaming


def test_export_transcriptome_respects_splicing_and_strand(tmp_path):
    """FASTA sequence is spliced in transcript orientation; GTF retains genomic coordinates."""
    genome=tmp_path/'genome.fa';genome.write_text('>chr1\nAAAACCCCGGGGTTTT\n');pysam.faidx(str(genome))
    ts=[Transcript('p','g','chr1','+',((0,4),(8,12))),Transcript('m','g','chr1','-',((0,4),(8,12)))]
    result=export_transcriptome(ts,genome,tmp_path/'models')
    text=Path(result['fasta']).read_text()
    assert '>p\nAAAAGGGG\n' in text
    assert '>m\nCCCCTTTT\n' in text
    assert '\texon\t1\t4\t' in Path(result['gtf']).read_text()
    with pytest.raises(FileExistsError):export_transcriptome(ts,genome,tmp_path/'models')


def test_subset_catalog_keeps_alias_structures_and_provenance(tmp_path):
    """Panel selection retains original PB identities and rejects ambiguous gene names."""
    data={'genes':{'g1':{'gene_id':'g1','gene_name':'A','chrom':'chr1','strand':'+','ref':[],'lr':[{'id':'PB.1','exons':[[0,20]],'fl':3}]},'g2':{'gene_id':'g2','gene_name':'B','chrom':'chr2','strand':'+','ref':[],'lr':[]}},'provenance':{'source':'input'}}
    path=tmp_path/'all.json.gz'
    with gzip.open(path,'wt') as f:json.dump(data,f)
    result=subset_catalog(path,['A'],tmp_path/'panel.json')
    panel=json.loads((tmp_path/'panel.json').read_text())
    assert list(panel['genes'])==['g1']
    assert panel['genes']['g1']['lr'][0]['id']=='PB.1'
    assert panel['provenance']['upstream']['source']=='input'
    assert result['gene_ids']==['g1']
    with pytest.raises(ValueError):subset_catalog(path,['missing'],tmp_path/'bad.json')


def test_short_read_panel_preserves_tags_and_forward_sequence(tmp_path):
    """Shared SR BAM and FASTQ contain the same eligible primary gene-assigned reads."""
    bam=tmp_path/'in.bam'
    with pysam.AlignmentFile(bam,'wb',header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000}]}) as out:
        for i,(cell,gene) in enumerate([('c','g'),('c','g'),('x','g'),('c','g;h')]):
            r=pysam.AlignedSegment();r.query_name=f'r{i}';r.query_sequence='AAAACCCC';r.query_qualities=pysam.qualitystring_to_array('IIIIIIII')
            r.flag=16;r.reference_id=0;r.reference_start=100+i;r.cigarstring='8M';r.mapping_quality=255
            r.set_tag('CB',cell);r.set_tag('UB',f'u{i}');r.set_tag('GX',gene);out.write(r)
    pysam.index(str(bam))
    t=Transcript('t','g','chr1','-',((90,200),))
    qc=extract_short_reads(bam,[t],{'c'},tmp_path/'panel')
    assert qc['accepted_reads']==2
    with pysam.AlignmentFile(tmp_path/'panel/sr.bam','rb') as f:
        reads=list(f);assert len(reads)==2 and reads[0].has_tag('GX')
    assert (tmp_path/'panel/pseudobulk.fastq').read_text().count('GGGGTTTT')==2
    assert Path(tmp_path/'panel/sr.bam.bai').exists()


def test_streaming_short_reads_match_window_selection_without_duplicate_fetches(tmp_path):
    """A genomic pass preserves window rules while overlapping loci write each record once."""
    bam=tmp_path/'in.bam'
    header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000},{'SN':'chr2','LN':500}]}
    rows=[
        ('keep_g',0,110,'g','c','u1'), ('keep_h',0,120,'h','c','u2'),
        ('antisense',16,130,'g','c','u4'),
        ('multi_gx',0,140,'g;h','c','u5'), ('bad_cell',0,150,'g','x','u6'),
        ('missing_umi',0,160,'g','c',None), ('wrong_locus',0,700,'g','c','u3'), ('unsupported_gene',0,800,'z','c','u7')]
    with pysam.AlignmentFile(bam,'wb',header=header) as out:
        for name,flag,start,gene,cell,umi in rows:
            r=pysam.AlignedSegment();r.query_name=name;r.query_sequence='AACCGGTT';r.query_qualities=pysam.qualitystring_to_array('IIIIIIII')
            r.flag=flag;r.reference_id=0 if start>=100 else 1;r.reference_start=start;r.cigarstring='8M';r.mapping_quality=255
            r.set_tag('CB',cell);r.set_tag('GX',gene)
            if umi is not None:r.set_tag('UB',umi)
            out.write(r)
    pysam.index(str(bam))
    transcripts=[Transcript('tg','g','chr1','+',((100,200),)),Transcript('th','h','chr1','+',((100,200),))]
    old=extract_short_reads(bam,transcripts,{'c'},tmp_path/'old')
    new=extract_short_reads_streaming(bam,transcripts,{'c'},tmp_path/'new')
    with pysam.AlignmentFile(tmp_path/'old/sr.bam','rb') as handle:old_names=[r.query_name for r in handle]
    with pysam.AlignmentFile(tmp_path/'new/sr.bam','rb') as handle:new_names=[r.query_name for r in handle]
    assert new_names==old_names==['keep_g','keep_h']
    assert new['accepted_reads']==old['accepted_reads']==2
    assert new['alignment_records']==len(rows)
    assert new['outside_gene_locus']==1
    assert new['unselected_gene']==1
    assert new['multi_gene_assignment']==1


def test_run_salmon_retains_duplicate_and_short_annotation_ids(tmp_path):
    """Real Salmon output retains exact duplicate models and explicit zero short models."""
    import csv,shutil
    from bulk2cell.preparation import run_salmon
    salmon=shutil.which('salmon') or ('bulk2cell/.tools/salmon/bin/salmon' if Path('bulk2cell/.tools/salmon/bin/salmon').exists() else None)
    if salmon is None:pytest.skip('Salmon integration binary unavailable')
    sequence='ACGTACGTTGCAAGTCGATCGTACGATGCTAGCATCGATCGTAGCTAGCATGCTAGCTAGCATCGATCGTAGCTACGATCGATGCTAGCATCGATC'
    fasta=tmp_path/'tx.fa';fasta.write_text(f'>unique\n{sequence}\n>duplicate\n{sequence}\n>short\nACGTACGTTGCAAGTCGATC\n')
    read=sequence[:50];fastqs=[]
    for part in range(2):
        fastq=tmp_path/f"reads{part}.fq";fastq.write_text("".join(f"@r{part}_{i}\n{read}\n+\n"+("I"*len(read))+"\n" for i in range(5)));fastqs.append(fastq)
    result=run_salmon(fasta,fastqs,tmp_path/'salmon',salmon,threads=2)
    with open(result['quant_sf']) as handle:rows={row['Name']:row for row in csv.DictReader(handle,delimiter='\t')}
    assert set(rows)=={'unique','duplicate','short'}
    assert rows['short']['TPM']=='0.000000' and rows['short']['NumReads']=='0.000'
