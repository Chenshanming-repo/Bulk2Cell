"""Check expanded fitting concurrency without changing scientific identities."""
import importlib.util
from pathlib import Path
import os
import time
import json
import pytest
from scipy.io import mmread
import numpy as np


def load():
    """Load the standalone execution extension."""
    path=Path(__file__).parents[1]/'benchmark/run_full_fit_workers.py'
    spec=importlib.util.spec_from_file_location('fit_workers',path)
    module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def worker_pid(index):
    """Keep tasks alive long enough to observe all requested processes."""
    time.sleep(.2)
    return os.getpid()


def noop():
    """Skip sample BAM initialization for the scheduler-only test."""


def test_twenty_processes_and_pending_bound(monkeypatch):
    """Use exactly20 children while retaining a40task submission bound."""
    module=load()
    monkeypatch.setattr(module.base.pipeline,'_initialize_worker',noop)
    sent=[];done=[]
    def jobs():
        """Track outstanding task count."""
        for i in range(80):
            sent.append(i)
            assert len(sent)-len(done)<=40
            yield i
    for result in module.completion_results(worker_pid,jobs(),20):done.append(result)
    assert len(set(done))==20
    with pytest.raises(ValueError):list(module.completion_results(worker_pid,[],21))


def test_real_bam_exact_science_and_resume(tmp_path):
    """Extra workers preserve outputs and reject inconsistent execution resumes."""
    from test_full_pipeline import fixture_args
    module=load()
    outputs=[]
    for name in ['original','expanded']:
        args=fixture_args(tmp_path,'quantify-full',tmp_path/name)
        command=['quantify-full','--catalog',args.catalog,'--bam',args.bam,'--barcodes',args.barcodes,'--out',args.out,'--workers','2','--tau','0','--window','50']
        if name=='original':assert module.base.main(['--']+command)==0
        else:
            assert module.main(['--fit-workers','20','--']+command)==0
            assert module.main(['--fit-workers','20','--']+command+['--resume'])==0
            with pytest.raises(ValueError,match='provenance'):
                module.main(['--fit-workers','19','--']+command+['--resume'])
        outputs.append(Path(args.out))
    for variant in module.base.pipeline.VARIANTS:
        a,b=[p/variant for p in outputs]
        for filename in ['isoforms.tsv','groups.tsv','membership.tsv','gene_qc.tsv','barcodes.tsv']:
            assert (a/filename).read_bytes()==(b/filename).read_bytes()
        for filename in ['isoform_counts.mtx.gz','group_counts.mtx.gz']:
            np.testing.assert_array_equal(mmread(a/filename).toarray(),mmread(b/filename).toarray())
        manifest=json.loads((b/'run.json').read_text())
        assert manifest['execution_controls']['workers']==20
        assert manifest['execution_controls']['fit_pending_limit']==40
