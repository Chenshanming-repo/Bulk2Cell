"""Orchestrate annotation, molecule likelihoods and auditable sparse output.

Input files are read-only. Results are written into a new directory and a
run.json completion record is emitted last; an interrupted directory must not
be interpreted as a completed run. Matrices are CELLS x FEATURES, unlike 10x.
"""
from collections import defaultdict, Counter
from pathlib import Path
import csv
import gzip
import hashlib
import importlib.metadata
import json
import resource
import time
import numpy as np
from scipy import sparse
from scipy.io import mmwrite
from . import adapters
from .inference import compatible, distance, fit_distance, quantify_gene
from . import __version__


def fingerprint(path):
    """Record size/mtime and hash small inputs; avoid a full scan of large BAMs."""
    p=Path(path).resolve(); stat=p.stat()
    record={'path':str(p),'bytes':stat.st_size,'mtime_ns':stat.st_mtime_ns}
    if stat.st_size <= 64*1024*1024:
        digest=hashlib.sha256()
        with p.open('rb') as f:
            for block in iter(lambda:f.read(1024*1024),b''): digest.update(block)
        record['sha256']=digest.hexdigest()
    else: record['hash_policy']='metadata_only_large_input'
    return record


def write_tsv(path, header, rows):
    """Write a UTF-8 tabular file with explicit column names and LF newlines."""
    with open(path,'w',newline='') as handle:
        writer=csv.writer(handle,delimiter='\t',lineterminator='\n')
        writer.writerow(header); writer.writerows(rows)


def _write_matrix(path, rows, cols, values, shape):
    """Serialize nonzero fractional expected counts in sparse Matrix Market form."""
    matrix=sparse.coo_matrix((values,(rows,cols)),shape=shape).tocsr()
    with gzip.open(path,'wb') as handle: mmwrite(handle,matrix,comment='cells x features; fractional expected UMI counts')


def run_quantification(args):
    """Execute a CLI namespace and return the JSON-serializable completion record."""
    started=time.monotonic()
    if Path(args.out).exists(): raise FileExistsError(f'output already exists: {args.out}')
    if args.catalog:
        transcripts,annotation_qc=adapters.load_regionquant_catalog(args.catalog)
    else:
        transcripts,annotation_qc=adapters.load_transcript_catalog(args.reference,args.isoseq_gff,args.classification)
    gene_aliases=annotation_qc.get('gene_aliases',{})
    selected={gene_aliases.get(g,g) for g in args.genes} if args.genes else {t.gene_id for t in transcripts}
    known={t.gene_id for t in transcripts}
    missing=selected-known
    if missing: raise ValueError('unknown gene IDs (use exact IDs): '+','.join(sorted(missing)))
    transcripts=[t for t in transcripts if t.gene_id in selected]
    if not transcripts: raise ValueError('no transcripts selected')
    barcodes=sorted(adapters.load_barcodes(args.barcodes))
    if not barcodes: raise ValueError('barcode whitelist is empty')
    aliases=annotation_qc.get('aliases',{})
    sr,sr_qc=adapters.load_salmon_quant(args.salmon,transcripts,aliases=aliases) if args.salmon else ({},{'fallback':'uniform_within_gene'})
    tes,tes_qc=adapters.load_tes_tsv(args.tes) if args.tes else ({},{'fallback':'annotation_TES'})
    canonical_tes=defaultdict(list)
    ids={t.id for t in transcripts}
    for key,points in tes.items():
        target=aliases.get(key,key)
        if target in ids: canonical_tes[target].extend(points)
    tes_qc['matched_transcripts']=len(canonical_tes)
    molecules,bam_qc=adapters.extract_molecules(args.bam,transcripts,selected,set(barcodes))
    by_gene=defaultdict(list); models=defaultdict(list)
    for t in transcripts: models[t.gene_id].append(t)
    for m in molecules: by_gene[m.gene_id].append(m)
    # Pool only structurally unique SR molecule distances. A bounded reservoir
    # avoids retaining all training values while preserving deterministic input.
    rng=np.random.default_rng(args.seed); training=[]; seen=0
    for gid,ts in models.items():
        if len(by_gene[gid])*len(ts)>args.max_likelihood_entries:
            raise ValueError(f'{gid}: likelihood budget exceeded; run smaller gene batches or increase --max-likelihood-entries')
        for m in by_gene[gid]:
            hits=[t for t in ts if compatible(m,t)]
            if len(hits)==1:
                d=distance(m,hits[0]); seen+=1
                if len(training)<2000: training.append(d)
                else:
                    j=int(rng.integers(seen))
                    if j<2000: training[j]=d
    model=fit_distance(training,args.bandwidth)
    out=Path(args.out); out.mkdir(parents=True)
    cell_index={c:i for i,c in enumerate(barcodes)}
    ir=[]; ic=[]; iv=[]; gr=[]; gc=[]; gv=[]
    features=[]; groups=[]; membership=[]; gene_qc=[]; totals=Counter()
    for gid in sorted(models):
        ts=sorted(models[gid],key=lambda t:t.id)
        result=quantify_gene(by_gene[gid],ts,sr,canonical_tes,tau=args.tau,window=args.window,
            model=model,em_strength=args.em_strength,max_iter=args.max_iter)
        ioffset=len(features); goffset=len(groups)
        for j,t in enumerate(ts): features.append((t.id,gid,t.source,t.lr_count,float(result['prior'][j])))
        for j,indices in enumerate(result['groups']):
            group_id=f'{gid}:group{j+1}'
            label='prior_informed_decomposition' if len(indices)>1 else 'single_member_group'
            groups.append((group_id,gid,len(indices),label))
            denom=result['prior'][indices].sum()
            for k in indices:
                weight=float(result['prior'][k]/denom) if denom else 1/len(indices)
                membership.append((group_id,ts[k].id,weight,label))
        for name,rr,cc,vv,offset in [('isoform_counts',ir,ic,iv,ioffset),('group_counts',gr,gc,gv,goffset)]:
            matrix=result[name]; r,c=np.nonzero(matrix)
            rr.extend(cell_index[result['cells'][i]] for i in r); cc.extend((c+offset).tolist()); vv.extend(matrix[r,c].tolist())
        qc=result['qc']; gene_qc.append((gid,len(ts),len(result['groups']),len(by_gene[gid]),qc['assigned_molecules'],qc['incompatible_molecules'],qc['em_nonconverged_cells'],qc['max_em_iterations']))
        for key in ('assigned_molecules','incompatible_molecules','em_nonconverged_cells'): totals[key]+=qc[key]
    if not np.isclose(sum(iv),totals['assigned_molecules']) or not np.isclose(sum(gv),totals['assigned_molecules']):
        raise RuntimeError('molecule conservation failed')
    _write_matrix(out/'isoform_counts.mtx.gz',ir,ic,iv,(len(barcodes),len(features)))
    _write_matrix(out/'group_counts.mtx.gz',gr,gc,gv,(len(barcodes),len(groups)))
    write_tsv(out/'barcodes.tsv',['barcode'],((c,) for c in barcodes))
    write_tsv(out/'isoforms.tsv',['transcript_id','gene_id','source','lr_count','hybrid_prior'],features)
    write_tsv(out/'groups.tsv',['group_id','gene_id','members','interpretation'],groups)
    write_tsv(out/'membership.tsv',['group_id','transcript_id','within_group_weight','interpretation'],membership)
    write_tsv(out/'gene_qc.tsv',['gene_id','isoforms','groups','input_molecules','assigned_molecules','incompatible_molecules','em_nonconverged_cells','max_em_iterations'],gene_qc)
    write_tsv(out/'distance_training.tsv',['distance'],((d,) for d in training))
    # Export the actual union for Salmon indexing and reproducible model review.
    with (out/'transcripts.gtf').open('w') as handle:
        for t in transcripts:
            for a,b in t.exons:
                handle.write(f'{t.chrom}\tbulk2cell\texon\t{a+1}\t{b}\t.\t{t.strand}\t.\tgene_id "{t.gene_id}"; transcript_id "{t.id}";\n')
    inputs=[getattr(args,k,None) for k in ('catalog','reference','isoseq_gff','classification','bam','barcodes','salmon','tes')]
    index=Path(args.bam+'.bai')
    if index.exists(): inputs.append(str(index))
    versions={p:importlib.metadata.version(p) for p in ('numpy','scipy','pysam')}
    manifest={'software':'bulk2cell','version':__version__,'status':'complete','matrix_orientation':'cells_by_features',
        'parameters':vars(args),'inputs':[fingerprint(p) for p in inputs if p],'dependencies':versions,
        'annotation_qc':annotation_qc,'salmon_qc':sr_qc,'tes_qc':tes_qc,'bam_qc':bam_qc,'qc':dict(totals),
        'cells':len(barcodes),'isoforms':len(features),'groups':len(groups),'genes':len(models),
        'distance_model':{'kind':'reflected_gaussian_kde' if training else 'exponential_fallback','bandwidth':args.bandwidth,'training_eligible_molecules':seen,'training_retained':len(training)},
        'elapsed_seconds':time.monotonic()-started,'peak_rss_mib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        'limitations':['GX single-gene assignments required; ambiguous or novel unassigned genes excluded',
            'terminal-window groups refined by observed likelihoods; not a universal identifiability proof',
            'within-group transcript estimates use bulk priors and cannot establish cell-specific switching',
            'per-cell EM uses prior regularization; group counts are expected, not uniquely observed counts']}
    with (out/'run.json').open('w') as handle: json.dump(manifest,handle,indent=2,allow_nan=False)
    return manifest
