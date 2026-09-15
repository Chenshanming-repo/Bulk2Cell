"""Compose preparation, FLNC alignment, and empirical TES command workflows."""
from __future__ import annotations
import json,struct
from pathlib import Path
from .adapters import load_barcodes,load_regionquant_catalog
from .flnc import extract_flnc_fasta
from .pipeline import fingerprint
from .preparation import export_transcriptome,extract_short_reads,run_logged,run_salmon,subset_catalog
from .tes import discover_isoseq_paths,estimate_tes,load_consensus_assignments,load_raw_assignments,write_tes_outputs

def _ids(transcripts,aliases):
    """Return retained canonical and original transcript identifiers."""
    canonical={x.id for x in transcripts};return canonical|{x for x,y in aliases.items() if y in canonical}

def estimate_tes_workflow(args):
    """Estimate separately labeled consensus or raw-read TES evidence."""
    tx,qc0=load_regionquant_catalog(args.catalog); aliases=qc0.get('aliases',{}); selected=_ids(tx,aliases)
    found=discover_isoseq_paths(args.isoseq_root) if args.isoseq_root else {}
    bam=Path(args.bam) if args.bam else (found.get('consensus_bam') if args.mode=='consensus' else None)
    table=Path(args.assignments) if args.assignments else found.get('group' if args.mode=='consensus' else 'read_stat')
    if bam is None or table is None: raise ValueError('mapped BAM and matching assignment table are required')
    loader=load_consensus_assignments if args.mode=='consensus' else load_raw_assignments
    assignments=loader(table,selected)
    rows,summary,qc=estimate_tes(bam,tx,assignments,mode=args.mode,aliases=aliases,selected_ids=selected,
        weight_mode=args.weight_mode,min_mapq=args.min_mapq,max_3prime_softclip=args.max_3prime_softclip,
        junction_tolerance=args.junction_tolerance,tes_window=args.tes_window)
    qc.update(catalog_qc=qc0,isoseq_root=str(args.isoseq_root) if args.isoseq_root else None,
        inputs={'catalog':fingerprint(args.catalog),'bam':fingerprint(bam),'assignments':fingerprint(table)})
    write_tes_outputs(args.out,rows,summary,qc);return qc

def extract_flnc_workflow(args):
    """Extract raw FLNC reads assigned to retained catalog transcripts."""
    out=Path(args.out); meta=Path(str(out)+'.run.json')
    if out.exists(): raise FileExistsError(str(out))
    if meta.exists() or (args.assignments_json and Path(args.assignments_json).exists()): raise FileExistsError(str(meta))
    tx,qc0=load_regionquant_catalog(args.catalog); selected=_ids(tx,qc0.get('aliases',{}))
    found=discover_isoseq_paths(args.isoseq_root) if args.isoseq_root else {}
    
    candidates=sorted(Path(args.isoseq_root).rglob('*.flnc.bam')) if args.isoseq_root else []
    if not args.bam and len(candidates)>1: raise ValueError(f'ambiguous raw FLNC BAM: {len(candidates)} matches')
    bam=Path(args.bam) if args.bam else (candidates[0] if candidates else None)
    table=Path(args.assignments) if args.assignments else found.get('read_stat')
    if bam is None or table is None: raise ValueError('raw FLNC BAM and read_stat assignments are required')
    assignments=load_raw_assignments(table,selected)
    qc=extract_flnc_fasta(bam,assignments,out,args.pbi,allow_missing=args.allow_missing)
    if args.assignments_json:
        path=Path(args.assignments_json)
        if path.exists(): raise FileExistsError(str(path))
        path.parent.mkdir(parents=True,exist_ok=True)
        path.write_text(json.dumps(assignments,sort_keys=True)+'\n')
    record=dict(qc,status='partial' if qc['missing'] else 'complete',evidence='raw_read',assignment_count=len(assignments),
        inputs={'catalog':fingerprint(args.catalog),'bam':fingerprint(bam),'assignments':fingerprint(table)})
    meta.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n');return record

def _validate_index(path):
    """Reject minimap indexes incompatible with pbmm2 ISOSEQ seeding."""
    path=Path(path)
    if path.suffix!='.mmi': return
    header=path.open('rb').read(24)
    if len(header)<24 or header[:4]!=b'MMI\x02': raise ValueError(f'unrecognized minimap2 index: {path}')
    w,k,_,flags=struct.unpack('<4i',header[4:20])
    if flags & 1: raise ValueError('ISOSEQ index must not use homopolymer-compressed minimizers')
    if flags & 2: raise ValueError('ISOSEQ index must retain reference sequences')
    if (w,k)!=(5,15): raise ValueError(f'index uses w={w},k={k}; ISOSEQ requires w=5,k=15; supply genome FASTA')

def align_flnc_workflow(args):
    """Align raw FLNC FASTA with pbmm2 ISOSEQ and record exact provenance."""
    out=Path(args.out);meta=Path(str(out)+'.run.json');log=Path(str(out)+'.log')
    if any(x.exists() for x in (out,meta,log,Path(str(log)+'.json'),Path(str(out)+'.bai'),Path(str(out)+'.csi'))): raise FileExistsError(str(out))
    if min(args.threads,args.sort_threads)<1: raise ValueError('thread counts must be positive')
    _validate_index(args.genome)
    cmd=[args.pbmm2,'align',args.genome,args.fasta,out,'--preset','ISOSEQ','--sort','--unmapped','-j',args.threads,'-J',args.sort_threads]
    run=run_logged(cmd,log)
    if not out.is_file(): raise RuntimeError(f'pbmm2 did not create {out}')
    # Exit status alone cannot establish that alignment output is consumable.
    import pysam
    try:
        pysam.quickcheck(str(out))
        with pysam.AlignmentFile(str(out),'rb') as bam:
            if not bam.has_index(): raise ValueError('sorted alignment lacks its BAM index')
            alignment_records=sum(1 for _ in bam.fetch(until_eof=True))
    except (OSError,ValueError,pysam.SamtoolsError) as error:
        raise ValueError(f'invalid pbmm2 BAM output: {error}') from error
    record=dict(run,status='complete',evidence='raw_read',preset='ISOSEQ',output_bam=str(out),alignment_records=alignment_records,
        inputs={'genome':fingerprint(args.genome),'fasta':fingerprint(args.fasta)})
    meta.write_text(json.dumps(record,indent=2,sort_keys=True)+'\n');return record

def prepare_workflow(args):
    """Create selected annotations, matched short reads, and optional Salmon priors."""
    out=Path(args.out)
    if args.threads<1: raise ValueError('threads must be positive')
    if out.exists(): raise FileExistsError(str(out))
    out.mkdir(); selection=subset_catalog(args.catalog,args.genes,out/'catalog.json')
    tx,qc=load_regionquant_catalog(out/'catalog.json');union=export_transcriptome(tx,args.genome,out/'union')
    reference=export_transcriptome([x for x in tx if 'reference' in x.source],args.genome,out/'reference')
    short=extract_short_reads(args.bam,tx,load_barcodes(args.barcodes),out/'short_reads',args.threads); salmon={}
    if args.salmon:
        for name in ('reference','union'): salmon[name]=run_salmon(out/f'{name}.fa',out/'short_reads/pseudobulk.fastq',out/f'{name}_salmon',args.salmon,args.threads)
    record={'status':'complete','selected_genes':selection['gene_ids'],'catalog_qc':qc,'short_read_qc':short,
        'transcriptomes':{'union':union,'reference':reference},'salmon':salmon,
        'inputs':{k:fingerprint(v) for k,v in {'catalog':args.catalog,'genome':args.genome,'bam':args.bam,'barcodes':args.barcodes}.items()}}
    (out/'run.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n');return record
