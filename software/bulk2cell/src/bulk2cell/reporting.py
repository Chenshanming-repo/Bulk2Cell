"""File-oriented clustering comparisons with explicit barcode alignment."""
from dataclasses import asdict
from pathlib import Path
import csv
import json
import numpy as np
from .evaluation import load_matrix_market, evaluate_cells, read_cell_ranger_clusters, compare_cell_matrices
from .pipeline import write_tsv


def read_column(path, column):
    """Read one named TSV column from bulk2cell metadata."""
    with open(path) as handle: return [r[column] for r in csv.DictReader(handle,delimiter='\t')]


def _finite_json(value):
    """Map undefined numerical metrics to JSON null rather than nonstandard NaN."""
    if isinstance(value,dict): return {k:_finite_json(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)): return [_finite_json(v) for v in value]
    if isinstance(value,(float,np.floating)) and not np.isfinite(value): return None
    return value


def evaluate_directory(args):
    """Evaluate a completed run, excluding zero-count cells and recording exclusions.

    Original SCALPEL matrices must be converted to the same documented file
    layout and cell-by-feature orientation. Feature sets may differ; retained
    barcode sets must match exactly, preventing positional label comparisons.
    """
    source=Path(args.run); out=Path(args.out)
    if out.exists(): raise FileExistsError(f'output already exists: {out}')
    if not (source/'run.json').exists(): raise ValueError('run directory lacks completion manifest')
    cells=read_column(source/'barcodes.tsv','barcode')
    features=read_column(source/'isoforms.tsv','transcript_id')
    matrix,_,_=load_matrix_market(source/'isoform_counts.mtx.gz',cells,features)
    labels=read_cell_ranger_clusters(args.clusters) if args.clusters else None
    keep=np.asarray(matrix.sum(axis=1)).ravel()>0
    active=[c for c,k in zip(cells,keep) if k]
    result=evaluate_cells(matrix[keep],active,args.n_clusters,reference_labels=labels,random_state=args.seed,compute_umap=args.umap)
    report={'cells_evaluated':len(active),'zero_count_cells_excluded':int((~keep).sum()),
        'method':result.method,'ari_cell_ranger_reference':result.adjusted_rand,
        'nmi_cell_ranger_reference':result.normalized_mutual_info,'silhouette':result.silhouette,
        'interpretation':'Cell Ranger agreement is not isoform accuracy ground truth',
        'parameters':vars(args)}
    if args.original:
        original=Path(args.original)
        bc=read_column(original/'barcodes.tsv','barcode'); ft=read_column(original/'isoforms.tsv','transcript_id')
        other,_,_=load_matrix_market(original/'isoform_counts.mtx.gz',bc,ft)
        index={c:i for i,c in enumerate(bc)}
        if any(c not in index for c in active): raise ValueError('original matrix missing evaluated barcodes')
        comparison=compare_cell_matrices(matrix[keep],active,other[[index[c] for c in active]],active,args.n_clusters,random_state=args.seed)
        report['original_scalpel_comparison']=asdict(comparison)
    out.mkdir(parents=True)
    write_tsv(out/'clusters.tsv',['barcode','cluster'],zip(active,result.clusters.tolist()))
    write_tsv(out/'embedding.tsv',['barcode']+[f'SVD{i+1}' for i in range(result.embedding.shape[1])],((c,*row) for c,row in zip(active,result.embedding.tolist())))
    if result.umap is not None: write_tsv(out/'umap.tsv',['barcode','UMAP1','UMAP2'],((c,*row) for c,row in zip(active,result.umap.tolist())))
    report=_finite_json(report)
    (out/'evaluation.json').write_text(json.dumps(report,indent=2,allow_nan=False))
    return report
