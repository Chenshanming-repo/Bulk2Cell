#!/usr/bin/env python3
"""Run frozen bulk2cell with bounded completion-order scheduling for fits only.

Usage: python run_full_completion_order.py [--adopt-existing
--transition-record FILE] -- quantify-full ... [--resume]

The original training scheduler, numerical functions, and identity checks are
unchanged. A separate, immutable OUT.scheduler.json binds this launcher and
its scheduling policy to the scientific identity. Resume of a pre-launcher
checkpoint requires explicit adoption plus a JSON transition record documenting
its verified snapshot and interrupted time. The caller performs that snapshot;
this launcher never imports checkpoints or rewrites their identities.

SIGINT/SIGTERM cancel pending fit tasks and terminate this launcher's workers.
Already committed SQLite gene results survive; in-flight and completed but
uncommitted results are recomputed on resume. SIGKILL cannot run this cleanup.
Use this same launcher for subsequent resumes so its provenance is enforced.
"""
from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from contextlib import contextmanager
import fcntl
import hashlib
import json
import multiprocessing
import os
from pathlib import Path
import signal
import sys
import tempfile

from bulk2cell import cli
from bulk2cell import full_pipeline as pipeline


def _worker_initialize():
    """Retain the original BAM/thread setup and let the parent handle interrupts."""
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    signal.signal(signal.SIGTERM, signal.SIG_DFL)
    pipeline._initialize_worker()


def _abort_pool(pool):
    """Terminate owned workers promptly on Python 3.10, then reap the executor.

    Python 3.10 has no public terminate_workers API. Snapshot its process table
    before shutdown; these are only children created by this executor.
    """
    processes = list(pool._processes.values())
    for process in processes:
        if process.is_alive():
            process.terminate()
    for process in processes:
        process.join(timeout=0.25)
    for process in processes:
        if process.is_alive():
            process.kill()
    pool.shutdown(wait=True, cancel_futures=True)
    for process in processes:
        process.join()


def completion_results(function, genes, workers):
    """Yield ready gene results with at most twice the workers outstanding."""
    if not 1 <= workers <= 8:
        raise ValueError('workers must be between 1 and 8')
    pool = ProcessPoolExecutor(max_workers=workers,
        mp_context=multiprocessing.get_context('fork'), initializer=_worker_initialize)
    finished = False
    try:
        iterator = iter(genes)
        pending = set()
        for _ in range(workers * 2):
            gid = next(iterator, None)
            if gid is not None:
                pending.add(pool.submit(function, gid))
        while pending:
            ready, _ = wait(pending, return_when=FIRST_COMPLETED)
            for future in ready:
                pending.remove(future)
                yield future.result()
                gid = next(iterator, None)
                if gid is not None:
                    pending.add(pool.submit(function, gid))
        finished = True
    finally:
        if finished:
            pool.shutdown(wait=True)
        else:
            _abort_pool(pool)


@contextmanager
def patched_scheduler():
    """Patch only fit scheduling, restoring the original dispatcher on exit."""
    original = pipeline._ordered_results
    def dispatch(function, genes, workers, pending_limit=None):
        """Keep training input order and its original pending-limit arguments."""
        if function is pipeline._fit_gene:
            yield from completion_results(function, genes, workers)
        else:
            yield from original(function, genes, workers, pending_limit=pending_limit)
    pipeline._ordered_results = dispatch
    try:
        yield
    finally:
        pipeline._ordered_results = original


def _sha256(path):
    """Hash small immutable source/provenance files without shell interpolation."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _transition(path):
    """Require a readable JSON transition record and bind its exact contents."""
    path = Path(path).resolve()
    if not isinstance(json.loads(path.read_text()), dict):
        raise ValueError('transition record must be a JSON object')
    return {'path': str(path), 'sha256': _sha256(path)}


def ensure_provenance(args, *, adopt_existing=False, transition_record=None):
    """Validate/create adjacent provenance without weakening production identity."""
    out = Path(args.out).resolve()
    path = Path(str(out) + '.scheduler.json')
    identity = pipeline._identity(args)
    if out.exists():
        if not args.resume:
            raise FileExistsError(f'output already exists: {out}; use --resume')
        saved = out/'identity.json'
        if not saved.exists() or json.loads(saved.read_text()) != identity:
            raise ValueError('original scientific identity differs; scheduler cannot adopt it')
    elif args.resume:
        raise ValueError('resume requires an existing output directory')
    if args.train_only:
        raise ValueError('this launcher is for fits; use the original CLI for train-only')
    previous = json.loads(path.read_text()) if path.exists() else None
    transition = None
    if previous and previous.get('transition'):
        transition = _transition(previous['transition']['path'])
        if transition != previous['transition']:
            raise ValueError('scheduler provenance transition record changed')
    if transition_record is not None:
        transition = _transition(transition_record)
    if out.exists() and previous is None and not (adopt_existing and transition is not None):
        raise ValueError('adopt existing checkpoints explicitly with --adopt-existing --transition-record')
    if adopt_existing and not out.exists():
        raise ValueError('adopt-existing requires an existing checkpoint directory')
    source = Path(pipeline.__file__).parent
    record = {'schema': 'bulk2cell.scheduler.v1', 'output': str(out),
        'launcher_path': str(Path(__file__).resolve()), 'launcher_sha256': _sha256(__file__),
        'module_sha256': {name: _sha256(source/name) for name in
            ('full_pipeline.py', 'inference.py', 'adapters.py', 'models.py', 'pipeline.py', 'cli.py')},
        'scientific_identity': identity, 'transition': transition,
        'policy': {'fit_order': 'completion', 'training_order': 'original',
                   'workers': args.workers, 'max_workers': 8, 'pending_limit': args.workers*2,
                   'max_pending': 16, 'thread_limit_per_worker': 1,
                   'interrupt': 'terminate_owned_workers; recompute_uncommitted_genes'}}
    if previous is not None:
        if previous != record:
            raise ValueError('scheduler provenance differs; refusing resume')
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir=path.parent, delete=False) as handle:
            temporary = handle.name
            json.dump(record, handle, indent=2, allow_nan=False)
            handle.write('\n')
        os.link(temporary, path)
    finally:
        if temporary is not None:
            os.unlink(temporary)
    return path


def _interrupt(signum, frame):
    """Convert parent termination into orderly Python cancellation and cleanup."""
    raise KeyboardInterrupt(f'received signal {signum}')


def main(argv=None):
    """Validate provenance, run the unchanged CLI, and attach execution provenance."""
    arguments = list(sys.argv[1:] if argv is None else argv)
    root = argparse.ArgumentParser(description=__doc__)
    root.add_argument('--adopt-existing', action='store_true')
    root.add_argument('--transition-record', type=Path)
    root.add_argument('command', nargs=argparse.REMAINDER)
    options = root.parse_args(arguments)
    command = options.command
    if command and command[0] == '--':
        command = command[1:]
    args = cli.parser().parse_args(command)
    if args.command != 'quantify-full':
        root.error('only quantify-full is supported')
    out = Path(args.out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    # Separate launcher lock protects sidecar creation; the original runner's
    # lock remains in force for all checkpoint and matrix mutations.
    with Path(str(out)+'.scheduler.lock').open('a') as handle:
        try:
            fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError('scheduler launcher already running for output') from error
        with (out.parent/('.'+out.name+'.lock')).open('a') as original_lock:
            try:
                fcntl.flock(original_lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
            except BlockingIOError as error:
                raise RuntimeError('original quantification still running for output') from error
        provenance = ensure_provenance(args, adopt_existing=options.adopt_existing,
                                       transition_record=options.transition_record)
        previous_signals = {sig: signal.signal(sig, _interrupt) for sig in (signal.SIGINT, signal.SIGTERM)}
        try:
            with patched_scheduler():
                status = cli.main(command)
            if status == 0:
                reference = {'path': str(provenance), 'sha256': _sha256(provenance)}
                for directory in [out] + [out/v for v in pipeline.VARIANTS]:
                    manifest = directory/'run.json'
                    record = json.loads(manifest.read_text())
                    record['scheduler_provenance'] = reference
                    record['execution_controls']['fit_result_order'] = 'completion'
                    pipeline._atomic_json(manifest, record)
            return status
        except KeyboardInterrupt:
            print('Interrupted: committed genes retained; uncommitted genes recompute on resume.', file=sys.stderr)
            return 130
        finally:
            for sig, handler in previous_signals.items():
                signal.signal(sig, handler)


if __name__ == '__main__':
    raise SystemExit(main())
