"""Command-line configuration and concise input-error reporting."""
import argparse
import json
import sys
from . import __version__


def parser():
    """Build the public CLI; scientific defaults are visible in --help."""
    root=argparse.ArgumentParser(prog='bulk2cell',description=__doc__)
    root.add_argument('--version',action='version',version=__version__)
    commands=root.add_subparsers(dest='command',required=True)
    q=commands.add_parser('quantify',help='Quantify an indexed Cell Ranger BAM',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    source=q.add_mutually_exclusive_group(required=True)
    source.add_argument('--catalog',help='Existing regionquant catalog JSON[.gz]')
    source.add_argument('--reference',help='Reference exon GTF[.gz]')
    q.add_argument('--isoseq-gff',help='Pigeon GTF-style GFF or GFF3')
    q.add_argument('--classification',help='Pigeon classification TSV for LR counts/QC')
    q.add_argument('--bam',required=True,help='Coordinate-sorted indexed BAM with CB/UB/GX')
    q.add_argument('--barcodes',required=True,help='Filtered barcodes TSV[.gz]')
    q.add_argument('--salmon',help='Salmon quant.sf against the matching union annotation')
    q.add_argument('--tes',help='LR endpoint TSV: transcript_id, tes, optional count')
    q.add_argument('--genes',nargs='+',help='Normalized gene IDs or unique gene names; default all catalog genes')
    q.add_argument('--out',required=True,help='New output directory')
    q.add_argument('--tau',type=float,default=10,help='SR prior equivalent LR count')
    q.add_argument('--window',type=int,default=600,help='Spliced 3-prime grouping window in nt')
    q.add_argument('--bandwidth',type=float,default=30,help='Gaussian distance bandwidth in nt')
    q.add_argument('--em-strength',type=float,default=1,help='Hybrid prior pseudocount strength in cell/group EM')
    q.add_argument('--max-iter',type=int,default=500,help='Maximum EM iterations per cell/gene')
    q.add_argument('--max-likelihood-entries',type=int,default=10000000,help='Per-gene dense likelihood allocation guard; 0 disables the limit')
    q.add_argument('--seed',type=int,default=0,help='Distance reservoir sampling seed')
    f=commands.add_parser('quantify-full',parents=[q],add_help=False,help='Bounded resumable all-gene quantification and ablations')
    f.add_argument('--workers',type=int,default=4,help='Concurrent genes (positive worker count)')
    f.add_argument('--train-only',action='store_true',help='Publish a global training cache before priors/TES are available')
    f.add_argument('--training-cache',help='Verified training.json from an identical evidence training run')
    f.add_argument('--resume',action='store_true',help='Resume only with identical input fingerprints and parameters')
    e=commands.add_parser('evaluate',help='Compare cell clustering with reference labels')
    e.add_argument('--run',required=True,help='Completed bulk2cell result directory')
    e.add_argument('--out',required=True,help='New evaluation directory')
    e.add_argument('--clusters',help='Cell Ranger Barcode/Cluster CSV')
    e.add_argument('--original',help='Original SCALPEL converted matrix directory')
    e.add_argument('--n-clusters',type=int,default=10)
    e.add_argument('--seed',type=int,default=0)
    e.add_argument('--umap',action='store_true',help='Also compute optional UMAP coordinates')
    t=commands.add_parser('estimate-tes',help='Estimate empirical TES from Iso-Seq alignments')
    t.add_argument('--catalog',required=True);t.add_argument('--isoseq-root');t.add_argument('--mode',choices=('consensus','raw'),default='consensus')
    t.add_argument('--bam');t.add_argument('--assignments');t.add_argument('--out',required=True)
    t.add_argument('--weight-mode',choices=('consensus','support'),default='consensus');t.add_argument('--min-mapq',type=int,default=20)
    t.add_argument('--max-3prime-softclip',type=int,default=20);t.add_argument('--junction-tolerance',type=int,default=10);t.add_argument('--tes-window',type=int,default=200)
    x=commands.add_parser('extract-flnc',help='Extract selected raw FLNC reads to FASTA')
    x.add_argument('--catalog',required=True);x.add_argument('--isoseq-root');x.add_argument('--bam');x.add_argument('--assignments');x.add_argument('--pbi')
    x.add_argument('--out',required=True);x.add_argument('--assignments-json');x.add_argument('--allow-missing',action='store_true')
    a=commands.add_parser('align-flnc',help='Align raw FLNC FASTA with pbmm2 ISOSEQ')
    a.add_argument('--fasta',required=True);a.add_argument('--genome',required=True);a.add_argument('--out',required=True);a.add_argument('--pbmm2',default='pbmm2')
    a.add_argument('--threads',type=int,default=4);a.add_argument('--sort-threads',type=int,default=2)
    p=commands.add_parser('prepare',help='Prepare matched annotation and short-read inputs')
    p.add_argument('--catalog',required=True);p.add_argument('--genes',nargs='+',required=True);p.add_argument('--genome',required=True)
    p.add_argument('--bam',required=True);p.add_argument('--barcodes',required=True);p.add_argument('--out',required=True);p.add_argument('--threads',type=int,default=2)
    p.add_argument('--salmon',help='Salmon executable; omit to skip Salmon')
    return root


def main(argv=None):
    """Execute a subcommand; return 2 for invalid input without hiding software bugs."""
    root=parser(); args=root.parse_args(argv)
    try:
        if args.command == 'evaluate':
            from .reporting import evaluate_directory
            print(json.dumps(evaluate_directory(args),indent=2)); return 0
        if args.command in {'estimate-tes','extract-flnc','align-flnc','prepare'}:
            from .orchestration import align_flnc_workflow,estimate_tes_workflow,extract_flnc_workflow,prepare_workflow
            workflows={'estimate-tes':estimate_tes_workflow,'extract-flnc':extract_flnc_workflow,'align-flnc':align_flnc_workflow,'prepare':prepare_workflow}
            print(json.dumps(workflows[args.command](args),indent=2)); return 0
        from .pipeline import run_quantification
        if bool(args.isoseq_gff) != bool(args.classification): raise ValueError('--isoseq-gff and --classification must be supplied together')
        if args.catalog and args.isoseq_gff: raise ValueError('--catalog cannot be combined with raw IsoSeq inputs')
        if args.max_likelihood_entries<0: raise ValueError('--max-likelihood-entries must be nonnegative (0 disables the limit)')
        if args.command == 'quantify-full':
            from .full_pipeline import run_full_quantification
            result=run_full_quantification(args)
        else:
            result=run_quantification(args)
        print(json.dumps({k:result[k] for k in ('status','cells','genes','isoforms','groups','qc','elapsed_seconds')},indent=2))
        return 0
    except (ValueError,OSError,ImportError,RuntimeError) as error:
        print(f'bulk2cell: {error}',file=sys.stderr); return 2
