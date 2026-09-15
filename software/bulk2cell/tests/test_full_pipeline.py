"""Bounded quantification must reproduce scalar output and reject stale resumes."""
import json
from pathlib import Path
import numpy as np
import pysam
import pytest
from scipy.io import mmread
from bulk2cell.cli import parser
from bulk2cell.pipeline import run_quantification


def fixture_args(tmp_path, command, out):
    """Create two expressed genes and one zero-evidence gene in a real BAM."""
    cp=tmp_path/'catalog.json'; bp=tmp_path/'cells.tsv'; bam=tmp_path/'reads.bam'
    if not cp.exists():
        genes={}
        for gid,start in [('z',100),('a',500),('empty',800)]:
            genes[gid]={'gene_id':gid,'gene_name':gid,'chrom':'chr1','strand':'+',
                'ref':[{'id':gid+'1','exons':[[start,start+150]]}],
                'lr':[{'id':gid+'2','exons':[[start+20,start+150]],'fl':3,'strict':True}]}
        cp.write_text(json.dumps({'genes':genes})); bp.write_text('cell-1\ncell-2\nempty-1\n')
        with pysam.AlignmentFile(bam,'wb',header={'HD':{'SO':'coordinate'},'SQ':[{'SN':'chr1','LN':1000}]}) as handle:
            for gid,start in [('z',100),('a',500)]:
                for i,offset in enumerate([0,90,95,100]):
                    r=pysam.AlignedSegment(); r.query_name=gid+str(i); r.query_sequence='A'*20
                    r.flag=0;r.reference_id=0;r.reference_start=start+offset;r.mapping_quality=255;r.cigarstring='20M'
                    r.set_tag('CB','cell-'+str(i%2+1));r.set_tag('UB',str(i));r.set_tag('GX',gid);handle.write(r)
        pysam.index(str(bam))
    return parser().parse_args([command,'--catalog',str(cp),'--bam',str(bam),'--barcodes',str(bp),'--out',str(out),'--window','50'])


def test_full_matches_scalar_and_resumes(tmp_path):
    """Global KDE, group definitions, empty features and counts agree exactly."""
    scalar=fixture_args(tmp_path,'quantify',tmp_path/'scalar'); run_quantification(scalar)
    full=fixture_args(tmp_path,'quantify-full',tmp_path/'full');full.workers=2
    from bulk2cell.full_pipeline import run_full_quantification
    result=run_full_quantification(full)
    for variant in ['empirical_tes','annotation_tes']:
        dest=Path(full.out)/variant
        for filename in ['isoforms.tsv','groups.tsv','membership.tsv','gene_qc.tsv','distance_training.tsv','barcodes.tsv']:
            assert (dest/filename).read_text()==(Path(scalar.out)/filename).read_text()
        for filename in ['isoform_counts.mtx.gz','group_counts.mtx.gz']:
            np.testing.assert_allclose(mmread(dest/filename).toarray(),mmread(Path(scalar.out)/filename).toarray(),atol=1e-12,rtol=1e-12)
    assert result['genes']==3
    with pytest.raises(FileExistsError): run_full_quantification(full)
    full.resume=True
    assert run_full_quantification(full)['status']=='complete'
    full.tau=11
    with pytest.raises(ValueError,match='fingerprint|parameters'): run_full_quantification(full)


def test_precomputed_likelihood_equivalence():
    """Reusing an annotation likelihood leaves every scalar inference result intact."""
    from bulk2cell.inference import quantify_gene,likelihoods,fit_distance
    from bulk2cell.models import Transcript,Molecule
    ts=[Transcript('a','g','chr1','+',((100,300),),0,'ref'),Transcript('b','g','chr1','+',((150,300),),4,'lr')]
    ms=[Molecule('c','u','g',((200,220),),()),Molecule('d','v','g',((110,130),),())]
    model=fit_distance([20,30]);matrix=likelihoods(ms,ts,model,{},normalize=True)
    expected=quantify_gene(ms,ts,{}, {},model=model)
    actual=quantify_gene(ms,ts,{}, {},model=model,precomputed_likelihood=matrix)
    for key in ['prior','group_counts','isoform_counts']:np.testing.assert_array_equal(actual[key],expected[key])
    assert actual['qc']==expected['qc'];assert actual['groups']==expected['groups']


def test_cached_likelihood_is_scalar_exact():
    """Repeated read shapes retain exact TES mixture, tails, strand, and incompatibility."""
    from bulk2cell import inference
    from bulk2cell.models import Transcript,Molecule
    from bulk2cell.full_pipeline import cached_likelihoods
    for strand in ['+','-']:
        ts=[Transcript('a','g','chr1',strand,((100,200),(300,450)),4),
            Transcript('b','g','chr1',strand,((100,220),(300,450)),0)]
        ms=[Molecule(str(c),str(i),'g',blocks,junctions) for c in range(3)
            for i,(blocks,junctions) in enumerate([(((310,330),),()),(((180,200),(300,320)),((200,300),)),(((50,70),),())])]
        tes={'a':[(t, w) for t,w in ([(440,2),(460,1)] if strand=='+' else [(90,2),(110,1)])]}
        for model in [inference.fit_distance([]),inference.fit_distance([1,10,80,9000])]:
            expected=inference.likelihoods(ms,ts,model,tes,normalize=True)
            np.testing.assert_array_equal(cached_likelihoods(ms,ts,model,tes),expected)


def test_missing_shard_is_recomputed_and_changed_input_rejected(tmp_path):
    """Completion markers never excuse missing genes, and modified inputs invalidate resume."""
    import sqlite3
    from bulk2cell.full_pipeline import run_full_quantification
    args=fixture_args(tmp_path,'quantify-full',tmp_path/'full');args.workers=1
    run_full_quantification(args)
    with sqlite3.connect(Path(args.out)/'checkpoints.sqlite') as db:
        db.execute("DELETE FROM fits WHERE gene='a'")
    args.resume=True
    assert run_full_quantification(args)['genes']==3
    with sqlite3.connect(Path(args.out)/'checkpoints.sqlite') as db:
        assert db.execute('SELECT COUNT(*) FROM fits').fetchone()[0]==3
    Path(args.barcodes).write_text('different-cell\n')
    with pytest.raises(ValueError,match='fingerprint'):run_full_quantification(args)


def test_batched_em_is_scalar_exact():
    """Independent EM stopping and zero-prior recovery agree across variable UMI counts."""
    from bulk2cell.inference import _em
    from bulk2cell.full_pipeline import batched_em
    rng=np.random.default_rng(47)
    for groups in [1,2,9]:
        matrices=[]
        for n in [1,1,1,3,3,10,10,2]:
            matrix=rng.uniform(0.001,1,(n,groups));matrix[matrix<0.3]=0
            matrix[:,0]+=1
            matrices.append(matrix)
        prior=rng.uniform(size=groups);prior[-1]=0
        if prior.sum():prior/=prior.sum()
        else:prior[:]=1
        for strength in [0,1,10]:
            expected=[_em(m,prior,strength,1e-7,50) for m in matrices]
            actual=batched_em(matrices,prior,strength,1e-7,50)
            for old,new in zip(expected,actual):
                np.testing.assert_array_equal(old[0],new[0]); assert old[1:]==new[1:]


def test_failed_resume_revokes_completion(tmp_path,monkeypatch):
    """A failed missing-shard recomputation cannot leave a stale complete manifest."""
    import sqlite3
    import bulk2cell.full_pipeline as module
    args=fixture_args(tmp_path,'quantify-full',tmp_path/'full');args.workers=1
    module.run_full_quantification(args)
    with sqlite3.connect(Path(args.out)/'checkpoints.sqlite') as db:db.execute("DELETE FROM fits WHERE gene='a'")
    def failure(gid):
        """Emulate an explicit worker failure for a requested gene."""
        raise RuntimeError('gene failed')
    monkeypatch.setattr(module,'_fit_gene',failure);args.resume=True
    with pytest.raises(RuntimeError,match='gene failed'):module.run_full_quantification(args)
    assert not (Path(args.out)/'run.json').exists()
    for variant in module.VARIANTS:assert not (Path(args.out)/variant/'run.json').exists()


def test_training_cache_is_independent_of_priors_and_bound_to_reads(tmp_path):
    """Training can finish before TES/Salmon and is reusable only for identical evidence."""
    from bulk2cell.full_pipeline import run_full_quantification
    args=fixture_args(tmp_path,'quantify-full',tmp_path/'training');args.workers=1
    args.train_only=True
    result=run_full_quantification(args)
    assert result['status']=='training_complete'
    fit=fixture_args(tmp_path,'quantify-full',tmp_path/'full');fit.workers=1
    fit.training_cache=str(Path(args.out)/'training.json')
    assert run_full_quantification(fit)['status']=='complete'
    changed=fixture_args(tmp_path,'quantify-full',tmp_path/'changed');changed.workers=1
    changed.training_cache=fit.training_cache;changed.seed=7
    with pytest.raises(ValueError,match='training cache'):run_full_quantification(changed)


def test_concurrent_runner_is_rejected(tmp_path):
    """A filesystem lock prevents two processes mutating one checkpoint database."""
    import fcntl
    from bulk2cell.full_pipeline import run_full_quantification
    args=fixture_args(tmp_path,'quantify-full',tmp_path/'full')
    with (tmp_path/'.full.lock').open('w') as handle:
        fcntl.flock(handle,fcntl.LOCK_EX|fcntl.LOCK_NB)
        with pytest.raises(RuntimeError,match='already running'):run_full_quantification(args)


def test_checkpoint_uses_nonblocking_wal(tmp_path):
    """Checkpoint readers can inspect progress while transactional WAL writes proceed."""
    import sqlite3
    from bulk2cell.full_pipeline import run_full_quantification
    args=fixture_args(tmp_path,'quantify-full',tmp_path/'full');args.workers=1
    run_full_quantification(args)
    with sqlite3.connect(Path(args.out)/'checkpoints.sqlite') as db:
        assert db.execute('PRAGMA journal_mode').fetchone()[0]=='wal'


def _variable_training_result(index):
    """Return deterministic distances after unequal delays to scramble completion order."""
    import time
    time.sleep(.03 if index % 16 == 0 else .0001)
    return index, list(range(index * 100, (index + 1) * 100))


def _initialize_empty_worker():
    """Allow scheduling tests to run real workers without opening a BAM."""


def test_larger_pending_window_preserves_ordered_reservoir(monkeypatch):
    """Submission depth cannot change the global reservoir or its random state."""
    import bulk2cell.full_pipeline as module
    monkeypatch.setattr(module,'_initialize_worker',_initialize_empty_worker)
    snapshots=[]
    for pending_limit in [16,128]:
        rng=np.random.default_rng(47);training=[];seen=0;order=[]
        for index,distances in module._ordered_results(_variable_training_result,range(160),4,pending_limit=pending_limit):
            order.append(index)
            for d in distances:
                seen+=1
                if len(training)<2000:training.append(d)
                else:
                    j=int(rng.integers(seen))
                    if j<2000:training[j]=d
        assert order==list(range(160))
        snapshots.append((seen,training,rng.bit_generator.state))
    assert snapshots[0]==snapshots[1]


def test_only_training_increases_pending_window(tmp_path,monkeypatch):
    """The full runner increases training lookahead while retaining the fit default."""
    import bulk2cell.full_pipeline as module
    original=module._ordered_results;calls=[]
    def observed(function,genes,workers,**kwargs):
        """Record execution controls while preserving real extraction and fitting."""
        calls.append((function.__name__,kwargs.get('pending_limit')))
        yield from original(function,genes,workers,**kwargs)
    monkeypatch.setattr(module,'_ordered_results',observed)
    args=fixture_args(tmp_path,'quantify-full',tmp_path/'full');args.workers=1
    module.run_full_quantification(args)
    assert calls==[('_train_gene',128),('_fit_gene',None)]
