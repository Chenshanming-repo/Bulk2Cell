"""Verify bounded completion scheduling, frozen identities, and exact fit outputs."""
import importlib.util
import json
import multiprocessing
from pathlib import Path
import time

import numpy as np
import pytest
from scipy.io import mmread

LAUNCHER = Path(__file__).parents[1] / 'benchmark/run_full_completion_order.py'


def load_launcher():
    """Load the standalone launcher without executing its CLI."""
    assert LAUNCHER.exists(), 'completion-order launcher must exist'
    spec = importlib.util.spec_from_file_location('completion_launcher', LAUNCHER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def noop():
    """Initialize scheduling-only workers without sample BAM access."""


def unequal(index):
    """Create an early straggler and inexpensive independent results."""
    time.sleep(0.5 if index == 0 else 0.002)
    return index


def sleeping(index):
    """Leave a running task available for interrupt cleanup verification."""
    time.sleep(30 if index else 0.01)
    return index


def test_completion_order_refills_without_exceeding_pending_budget(monkeypatch):
    """Completed later jobs unlock new work while the initial slow job runs."""
    launcher = load_launcher()
    import bulk2cell.full_pipeline as pipeline
    monkeypatch.setattr(pipeline, '_initialize_worker', noop)
    consumed = []
    yielded = []
    def genes():
        """Observe submitted-but-not-consumed task count from the parent."""
        for index in range(40):
            consumed.append(index)
            assert len(consumed) - len(yielded) <= 4
            yield index
    for result in launcher.completion_results(unequal, genes(), 2):
        yielded.append(result)
    assert sorted(yielded) == list(range(40))
    assert yielded.index(0) > 4


def test_only_fit_dispatch_changes_and_training_arguments_survive(monkeypatch):
    """Training continues through its original ordered reservoir scheduler."""
    launcher = load_launcher()
    import bulk2cell.full_pipeline as pipeline
    calls = []
    def original(function, genes, workers, pending_limit=None):
        """Record original dispatch without mutating numerical inputs."""
        calls.append((function, workers, pending_limit))
        yield from genes
    monkeypatch.setattr(pipeline, '_ordered_results', original)
    with launcher.patched_scheduler():
        assert list(pipeline._ordered_results(pipeline._train_gene, [1, 2], 3, pending_limit=128)) == [1, 2]
    assert pipeline._ordered_results is original
    assert calls == [(pipeline._train_gene, 3, 128)]


def test_close_terminates_running_workers_promptly(monkeypatch):
    """Closing an interrupted scheduler cannot wait for long in-flight genes."""
    launcher = load_launcher()
    import bulk2cell.full_pipeline as pipeline
    monkeypatch.setattr(pipeline, '_initialize_worker', noop)
    before = {p.pid for p in multiprocessing.active_children()}
    results = launcher.completion_results(sleeping, range(4), 2)
    assert next(results) == 0
    started = time.monotonic()
    results.close()
    assert time.monotonic() - started < 5
    assert {p.pid for p in multiprocessing.active_children()} <= before


def test_provenance_precedes_output_and_rejects_changed_launcher_record(tmp_path):
    """Adjacent provenance is enforced without creating the output directory."""
    launcher = load_launcher()
    from test_full_pipeline import fixture_args
    args = fixture_args(tmp_path, 'quantify-full', tmp_path/'out')
    path = launcher.ensure_provenance(args)
    assert not Path(args.out).exists()
    record = json.loads(path.read_text())
    assert record['policy']['fit_order'] == 'completion'
    assert record['scientific_identity']['code']['full_pipeline.py']
    assert record['launcher_sha256']
    record['launcher_sha256'] = 'modified'
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match='provenance'):
        launcher.ensure_provenance(args)


def test_resume_requires_explicit_adoption_and_preserves_original_identity(tmp_path):
    """Existing runs need an explicit transition record; identity bypass is forbidden."""
    launcher = load_launcher()
    from test_full_pipeline import fixture_args
    import bulk2cell.full_pipeline as pipeline
    args = fixture_args(tmp_path, 'quantify-full', tmp_path/'old')
    Path(args.out).mkdir()
    identity = Path(args.out)/'identity.json'
    identity.write_text(json.dumps(pipeline._identity(args)))
    original = identity.read_bytes()
    args.resume = True
    with pytest.raises(ValueError, match='adopt'):
        launcher.ensure_provenance(args)
    transition = tmp_path/'transition.json'
    transition.write_text(json.dumps({'original_run': str(args.out), 'reason': 'resource-only scheduler'}))
    launcher.ensure_provenance(args, adopt_existing=True, transition_record=transition)
    assert identity.read_bytes() == original
    args.tau += 1
    with pytest.raises(ValueError, match='identity'):
        launcher.ensure_provenance(args)


def test_all_variants_are_scientifically_exact(tmp_path):
    """Ordered and completion-order fits publish identical counts, memberships and QC."""
    launcher = load_launcher()
    from test_full_pipeline import fixture_args
    import bulk2cell.full_pipeline as pipeline
    ordered = fixture_args(tmp_path, 'quantify-full', tmp_path/'ordered')
    reordered = fixture_args(tmp_path, 'quantify-full', tmp_path/'reordered')
    ordered.workers = reordered.workers = 2
    pipeline.run_full_quantification(ordered)
    with launcher.patched_scheduler():
        pipeline.run_full_quantification(reordered)
    for variant in pipeline.VARIANTS:
        left, right = Path(ordered.out)/variant, Path(reordered.out)/variant
        for filename in ['isoforms.tsv', 'groups.tsv', 'membership.tsv', 'gene_qc.tsv',
                         'distance_training.tsv', 'barcodes.tsv', 'transcripts.gtf']:
            assert (left/filename).read_bytes() == (right/filename).read_bytes()
        for filename in ['isoform_counts.mtx.gz', 'group_counts.mtx.gz']:
            a, b = mmread(left/filename).tocsr(), mmread(right/filename).tocsr()
            np.testing.assert_array_equal(a.indptr, b.indptr)
            np.testing.assert_array_equal(a.indices, b.indices)
            np.testing.assert_array_equal(a.data, b.data)
        a, b = [json.loads((p/'run.json').read_text()) for p in (left, right)]
        for key in ['qc', 'bam_qc', 'distance_model', 'variant_configuration', 'annotation_qc', 'tes_qc', 'salmon_qc']:
            assert a[key] == b[key]


def test_cli_attaches_enforced_scheduler_provenance(tmp_path):
    """The complete CLI records the immutable sidecar and retains ordinary resume checks."""
    launcher = load_launcher()
    from test_full_pipeline import fixture_args
    args = fixture_args(tmp_path, 'quantify-full', tmp_path/'cli')
    command = ['--', 'quantify-full', '--catalog', args.catalog, '--bam', args.bam,
               '--barcodes', args.barcodes, '--out', args.out, '--workers', '2', '--window', '50']
    assert launcher.main(command) == 0
    sidecar = Path(str(Path(args.out).resolve())+'.scheduler.json')
    for variant in ['', 'empirical_tes', 'annotation_tes', 'sr_only']:
        record = json.loads((Path(args.out)/variant/'run.json').read_text())
        assert record['scheduler_provenance']['path'] == str(sidecar)
        assert record['execution_controls']['fit_result_order'] == 'completion'
    assert launcher.main(command + ['--resume']) == 0
    with pytest.raises(ValueError, match='identity'):
        launcher.main(command + ['--resume', '--tau', '11'])


def test_sigterm_cleans_owned_fit_workers_and_preserves_checkpoint(tmp_path):
    """Actual SIGTERM exits promptly with committed training retained and no fit completion."""
    import os
    import subprocess
    import sys
    from test_full_pipeline import fixture_args
    args = fixture_args(tmp_path, 'quantify-full', tmp_path/'interrupted')
    program = tmp_path/'signal_fixture.py'
    program.write_text('''import importlib.util, os, sys, time
from pathlib import Path
spec=importlib.util.spec_from_file_location('launcher', sys.argv[1])
launcher=importlib.util.module_from_spec(spec);spec.loader.exec_module(launcher)
def slow_fit(gid):
    Path(sys.argv[2]+str(os.getpid())).write_text(gid)
    time.sleep(60)
launcher.pipeline._fit_gene=slow_fit
raise SystemExit(launcher.main(sys.argv[3:]))
''')
    marker = str(tmp_path/'worker-')
    command = [sys.executable, str(program), str(LAUNCHER), marker, '--', 'quantify-full',
               '--catalog', args.catalog, '--bam', args.bam, '--barcodes', args.barcodes,
               '--out', args.out, '--workers', '2', '--window', '50']
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        deadline = time.monotonic()+10
        while not list(tmp_path.glob('worker-*')) and process.poll() is None and time.monotonic()<deadline:
            time.sleep(0.02)
        markers = list(tmp_path.glob('worker-*'))
        assert markers, 'worker did not reach fit stage'
        started = time.monotonic()
        process.terminate()
        stdout, stderr = process.communicate(timeout=5)
        assert process.returncode == 130, (stdout, stderr)
        assert time.monotonic()-started < 5
        assert not (Path(args.out)/'run.json').exists()
        import sqlite3
        with sqlite3.connect(Path(args.out)/'checkpoints.sqlite') as db:
            assert db.execute("SELECT COUNT(*) FROM state WHERE key='training'").fetchone()[0] == 1
        for path in tmp_path.glob('worker-*'):
            with pytest.raises(ProcessLookupError):
                os.kill(int(path.name.split('-')[-1]), 0)
    finally:
        if process.poll() is None:
            process.kill()
            process.wait()
