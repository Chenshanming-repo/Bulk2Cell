"""Prepare matched annotations, short reads and external-tool inputs.

All source data are read-only. New output files/directories are required. FASTA
sequences follow transcript orientation; GTF intervals remain genomic and
one-based closed. Pseudobulk FASTQ contains exactly the accepted BAM records,
in their original sequencing orientation, without barcode/UMI sequence.
"""
from concurrent.futures import ProcessPoolExecutor
from collections import Counter, defaultdict
from pathlib import Path
import gzip
import json
import subprocess
import time
import pysam
from .adapters import normalize_gene_id
from .pipeline import fingerprint


def subset_catalog(path, genes, output_path):
    """Select exact gene IDs or unique symbols, preserving PB IDs and provenance."""
    output=Path(output_path)
    if output.exists(): raise FileExistsError(str(output))
    opener=gzip.open if str(path).endswith('.gz') else open
    with opener(path,'rt') as handle: payload=json.load(handle)
    available=payload.get('genes',{})
    if not isinstance(available,dict) or not available: raise ValueError('catalog must contain genes')
    selected=[]
    for token in genes:
        hits=[token] if token in available else [k for k,v in available.items() if v.get('gene_name')==token]
        if len(hits)!=1: raise ValueError(f'gene is missing or ambiguous: {token}')
        if hits[0] not in selected: selected.append(hits[0])
    if not selected: raise ValueError('at least one gene must be selected')
    result={'schema':payload.get('schema','regionquant.catalog.v1'),
        'genes':{gid:available[gid] for gid in selected},
        'provenance':{'upstream':payload.get('provenance',{}),'catalog_input':fingerprint(path),'selected_genes':selected}}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2))
    return {'catalog':str(output),'gene_ids':selected}


def export_transcriptome(transcripts, genome_path, output_prefix):
    """Export GTF and spliced FASTA against an existing indexed reference genome."""
    prefix=Path(output_prefix); gtf=Path(str(prefix)+'.gtf'); fasta=Path(str(prefix)+'.fa')
    if gtf.exists() or fasta.exists(): raise FileExistsError(str(prefix))
    if not Path(str(genome_path)+'.fai').exists(): raise ValueError('genome FASTA must already have a .fai index')
    if not transcripts or len({t.id for t in transcripts})!=len(transcripts): raise ValueError('transcripts must have unique IDs and be nonempty')
    for t in transcripts:
        if any(c in t.id+t.gene_id for c in '\t\n\r" '): raise ValueError('transcript/gene IDs cannot contain whitespace or quotes')
    prefix.parent.mkdir(parents=True,exist_ok=True)
    complement=str.maketrans('ACGTRYMKBDHVNacgtrymkbdhvn','TGCAYRKMVHDBNtgcayrkmvhdbn')
    with pysam.FastaFile(str(genome_path)) as genome,gtf.open('w') as annotation,fasta.open('w') as sequences:
        for t in transcripts:
            if t.chrom not in genome.references or t.exons[-1][1]>genome.get_reference_length(t.chrom):
                raise ValueError(f'transcript {t.id} exceeds reference contig bounds')
            sequence=''.join(genome.fetch(t.chrom,a,b) for a,b in t.exons).upper()
            if t.strand=='-': sequence=sequence.translate(complement)[::-1]
            sequences.write(f'>{t.id}\n{sequence}\n')
            attributes=f'gene_id "{t.gene_id}"; transcript_id "{t.id}"; gene_name "{t.gene_id}"; transcript_name "{t.id}";'
            for a,b in t.exons:
                annotation.write(f'{t.chrom}\tbulk2cell\texon\t{a+1}\t{b}\t.\t{t.strand}\t.\t{attributes}\n')
    return {'gtf':str(gtf),'fasta':str(fasta),'transcripts':len(transcripts)}


def extract_short_reads(bam_path, transcripts, barcodes, output_dir, threads=2):
    """Write a shared SR panel BAM and pseudobulk FASTQ with explicit selection QC.

    Corrected CB/UB tags, one GX, primary sense alignment and a selected gene
    locus are required. Overlapping fetch windows cannot duplicate BAM records.
    PCR records are retained here so each method can apply its own declared
    molecule handling; Salmon sees the same RNA sequence records.
    """
    out=Path(output_dir)
    if out.exists(): raise FileExistsError(str(out))
    by_gene=defaultdict(list)
    for t in transcripts: by_gene[normalize_gene_id(t.gene_id)].append(t)
    if not by_gene or not barcodes: raise ValueError('transcripts and barcode whitelist must be nonempty')
    loci={}
    for gid,models in by_gene.items():
        places={(t.chrom,t.strand) for t in models}
        if len(places)!=1: raise ValueError(f'gene {gid} spans loci')
        chrom,strand=next(iter(places));loci[gid]=(chrom,strand,min(t.exons[0][0] for t in models),max(t.exons[-1][1] for t in models))
    out.mkdir(parents=True)
    qc=Counter();seen=set();unsorted=out/'unsorted.bam'
    with pysam.AlignmentFile(str(bam_path),'rb') as bam,pysam.AlignmentFile(str(unsorted),'wb',header=bam.header) as writer,(out/'pseudobulk.fastq').open('w') as fastq:
        if not bam.has_index(): raise ValueError('source SR BAM must be indexed')
        for gid,(chrom,strand,start,end) in sorted(loci.items()):
            for read in bam.fetch(chrom,start,end):
                qc['fetch_events']+=1
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
                    qc['nonprimary_or_qcfail']+=1;continue
                if not read.has_tag('CB') or read.get_tag('CB') not in barcodes:
                    qc['barcode_excluded']+=1;continue
                if not read.has_tag('UB') or not read.get_tag('UB'):
                    qc['missing_umi']+=1;continue
                assigned=[normalize_gene_id(g) for g in str(read.get_tag('GX')).split(';') if g] if read.has_tag('GX') else []
                if assigned!=[gid]: qc['gene_assignment_excluded']+=1;continue
                if read.is_reverse!=(strand=='-'):qc['antisense']+=1;continue
                key=(read.query_name,read.flag,read.reference_id,read.reference_start,read.cigarstring)
                if key in seen:qc['repeated_fetch']+=1;continue
                seen.add(key)
                sequence=read.get_forward_sequence()
                if not sequence:qc['missing_sequence']+=1;continue
                qualities=read.get_forward_qualities()
                quality=pysam.array_to_qualitystring(qualities) if qualities is not None else 'I'*len(sequence)
                writer.write(read)
                fastq.write(f'@{read.query_name}\n{sequence}\n+\n{quality}\n')
                qc['accepted_reads']+=1
    pysam.sort('-@',str(threads),'-o',str(out/'sr.bam'),str(unsorted));pysam.index(str(out/'sr.bam'));unsorted.unlink()
    (out/'barcodes.tsv').write_text(''.join(c+'\n' for c in sorted(barcodes)))
    (out/'selection.json').write_text(json.dumps(dict(qc),indent=2))
    return dict(qc)


def run_logged(command, log_path):
    """Run an argument vector without a shell and record timing, exit code and log."""
    log=Path(log_path);log.parent.mkdir(parents=True,exist_ok=True)
    started=time.monotonic()
    with log.open('w') as handle:
        completed=subprocess.run([str(x) for x in command],stdout=handle,stderr=subprocess.STDOUT,check=False)
    record={'command':[str(x) for x in command],'elapsed_seconds':time.monotonic()-started,'returncode':completed.returncode,'log':str(log)}
    Path(str(log)+'.json').write_text(json.dumps(record,indent=2))
    if completed.returncode: raise RuntimeError(f'command failed ({completed.returncode}); see {log}')
    return record


def run_salmon(fasta, fastq, output_dir, executable='salmon', threads=4):
    """Build a transcriptome index and quantify unstranded single-end pseudobulk.

    The input is the sequenced RNA mate extracted from a 10x genome BAM, not the
    barcode/UMI mate. Fragment length defaults and Salmon regularization are
    explicit. A transcriptome-only targeted index is a functional benchmark;
    decoy-aware whole-transcriptome mapping is needed for production assessment.
    """
    out=Path(output_dir)
    if out.exists():raise FileExistsError(str(out))
    out.mkdir(parents=True)
    fastqs=[fastq] if isinstance(fastq,(str,Path)) else list(fastq)
    if not fastqs:raise ValueError("at least one FASTQ is required")
    commands=[run_logged([executable,'index','-t',fasta,'-i',out/'index','--keepDuplicates','-p',threads],out/'index.log'),
        run_logged([executable,'quant','-i',out/'index','-l','U','-r',*fastqs,'--fldMean','200','--fldSD','80',
            '--minScoreFraction','0.7','--vbPrior','5','--perNucleotidePrior','-p',threads,'-o',out/'quant'],out/'quant.log')]
    with pysam.FastxFile(str(fasta)) as handle: annotation_ids={entry.name for entry in handle}
    with (out/'quant/quant.sf').open() as handle:
        header=handle.readline().rstrip("\n").split("\t"); name_column=header.index("Name")
        quantified_ids={line.rstrip("\n").split("\t")[name_column] for line in handle if line.strip()}
    missing=sorted(annotation_ids-quantified_ids); unexpected=sorted(quantified_ids-annotation_ids)
    if missing or unexpected: raise RuntimeError(f'Salmon transcript ID coverage mismatch: missing={len(missing)}, unexpected={len(unexpected)}')
    record={'commands':commands,'fasta':fingerprint(fasta),'fastq':fingerprint(fastqs[0]) if len(fastqs)==1 else None,'fastqs':[fingerprint(path) for path in fastqs],
        'quant_sf':str(out/'quant/quant.sf'),'library_type':'U','fragment_mean':200,'fragment_sd':80,
        'keep_duplicates':True,'annotation_transcript_ids':len(annotation_ids),'quantified_transcript_ids':len(quantified_ids),
        'missing_transcript_ids':missing,'unexpected_transcript_ids':unexpected}
    (out/'run.json').write_text(json.dumps(record,indent=2))
    return record


def _extract_short_read_contig(task):
    """Extract one coordinate-sorted contig shard for bounded parallel preparation."""
    bam_path,chrom,loci,barcodes,bam_output,fastq_output=task
    qc=Counter()
    with pysam.AlignmentFile(str(bam_path),'rb') as bam, \
            pysam.AlignmentFile(str(bam_output),'wb',header=bam.header) as writer, \
            Path(fastq_output).open('w') as fastq:
        for read in bam.fetch(chrom):
            qc['alignment_records']+=1
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
                qc['nonprimary_or_qcfail']+=1;continue
            assigned=[normalize_gene_id(g) for g in str(read.get_tag('GX')).split(';') if g] if read.has_tag('GX') else []
            if len(assigned)!=1:
                qc['multi_gene_assignment' if len(assigned)>1 else 'missing_gene_assignment']+=1;continue
            locus=loci.get(assigned[0])
            if locus is None:
                qc['unselected_gene']+=1;continue
            locus_chrom,strand,start,end=locus
            if locus_chrom!=chrom or read.reference_end is None or read.reference_start>=end or read.reference_end<=start:
                qc['outside_gene_locus']+=1;continue
            if not read.has_tag('CB') or read.get_tag('CB') not in barcodes:
                qc['barcode_excluded']+=1;continue
            if not read.has_tag('UB') or not read.get_tag('UB'):
                qc['missing_umi']+=1;continue
            if read.is_reverse!=(strand=='-'):
                qc['antisense']+=1;continue
            sequence=read.get_forward_sequence()
            if not sequence:
                qc['missing_sequence']+=1;continue
            qualities=read.get_forward_qualities()
            quality=pysam.array_to_qualitystring(qualities) if qualities is not None else 'I'*len(sequence)
            writer.write(read);fastq.write(f'@{read.query_name}\n{sequence}\n+\n{quality}\n');qc['accepted_reads']+=1
    return dict(qc)


def extract_short_reads_streaming(bam_path, transcripts, barcodes, output_dir, threads=2):
    """Extract selected reads in bounded per-contig passes without global read identities.

    Each source contig is scanned once and yields a coordinate-sorted shard. Shards
    are concatenated in BAM header order, so no memory-heavy global set or resort is
    needed. Selection matches :func:`extract_short_reads`: primary sense records
    overlapping the single assigned GX gene window with whitelisted CB and nonempty
    UB. Unsupported catalog contigs are retained in explicit QC.
    """
    out=Path(output_dir)
    if out.exists():raise FileExistsError(str(out))
    if not transcripts or not barcodes:raise ValueError('transcripts and barcode whitelist must be nonempty')
    by_gene=defaultdict(list)
    for transcript in transcripts:by_gene[normalize_gene_id(transcript.gene_id)].append(transcript)
    loci={}
    for gid,models in by_gene.items():
        places={(t.chrom,t.strand) for t in models}
        if len(places)!=1:raise ValueError(f'gene {gid} spans loci')
        chrom,strand=next(iter(places));loci[gid]=(chrom,strand,min(t.exons[0][0] for t in models),max(t.exons[-1][1] for t in models))
    with pysam.AlignmentFile(str(bam_path),'rb') as bam:
        if not bam.has_index():raise ValueError('source SR BAM must be indexed')
        references=list(bam.references)
    unsupported=sorted(gid for gid,(chrom,_,_,_) in loci.items() if chrom not in set(references))
    supported_loci={gid:locus for gid,locus in loci.items() if gid not in set(unsupported)}
    used_contigs=set(locus[0] for locus in supported_loci.values())
    contigs=[chrom for chrom in references if chrom in used_contigs]
    out.mkdir(parents=True);shards=out/'shards';shards.mkdir()
    tasks=[(str(bam_path),chrom,supported_loci,set(barcodes),str(shards/f'{index:04d}.bam'),str(shards/f'{index:04d}.fastq')) for index,chrom in enumerate(contigs)]
    worker_count=max(1,min(int(threads),len(tasks)))
    if worker_count==1:results=[_extract_short_read_contig(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=worker_count) as pool:results=list(pool.map(_extract_short_read_contig,tasks))
    qc=Counter()
    for result in results:qc.update(result)
    bam_shards=[str(shards/f'{index:04d}.bam') for index in range(len(tasks))]
    if not bam_shards:raise ValueError('no catalog genes use source BAM contigs')
    pysam.cat('-o',str(out/'sr.bam'),*bam_shards);pysam.index(str(out/'sr.bam'))
    with (out/'pseudobulk.fastq').open('wb') as target:
        for index in range(len(tasks)):
            with (shards/f'{index:04d}.fastq').open('rb') as source:
                while chunk:=source.read(1024*1024):target.write(chunk)
    for path in shards.iterdir():path.unlink()
    shards.rmdir()
    (out/'barcodes.tsv').write_text(''.join(c+'\n' for c in sorted(barcodes)))
    report=dict(qc);report.update({'workers':worker_count,'source_contigs_scanned':contigs,'catalog_genes':len(loci),'supported_genes':len(supported_loci),'unsupported_contig_genes':unsupported})
    (out/'selection.json').write_text(json.dumps(report,indent=2))
    return report
