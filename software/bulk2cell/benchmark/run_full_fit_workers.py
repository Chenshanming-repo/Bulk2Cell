#!/usr/bin/env python3
"""Extend the frozen fitting launcher to at most20 single-thread fit workers.

The inner --workers remains its legacy admission/training control (at most8).
--fit-workers controls the actual fit pool only. Scientific identities and old
scheduler sidecars remain intact; OUT.fitworkers.json records effective controls.
Use this extension for all subsequent resumes of an extended run.
"""
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from contextlib import contextmanager
import argparse
import fcntl
import json
import multiprocessing
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parent))
import run_full_completion_order as base


def completion_results(function, genes, workers):
    """Yield completed fits with at most20 processes and twice as many tasks."""
    if not 1 <= workers <= 20:
        raise ValueError('fit workers must be between1 and20')
    pool=ProcessPoolExecutor(max_workers=workers,
        mp_context=multiprocessing.get_context('fork'),initializer=base._worker_initialize)
    finished=False
    try:
        iterator=iter(genes)
        pending=set()
        for _ in range(workers*2):
            gid=next(iterator,None)
            if gid is not None:pending.add(pool.submit(function,gid))
        while pending:
            ready,_=wait(pending,return_when=FIRST_COMPLETED)
            for future in ready:
                pending.remove(future)
                yield future.result()
                gid=next(iterator,None)
                if gid is not None:pending.add(pool.submit(function,gid))
        finished=True
    finally:
        if finished:pool.shutdown(wait=True)
        else:base._abort_pool(pool)


@contextmanager
def expanded_scheduler(workers):
    """Override only fit dispatch; retain original training and restore state."""
    original=base.pipeline._ordered_results
    def dispatch(function,genes,legacy_workers,pending_limit=None):
        """Select effective fit workers without altering training controls."""
        if function is base.pipeline._fit_gene:
            yield from completion_results(function,genes,workers)
        else:
            yield from original(function,genes,legacy_workers,pending_limit=pending_limit)
    base.pipeline._ordered_results=dispatch
    try:yield
    finally:base.pipeline._ordered_results=original


def main(argv=None):
    """Bind effective execution provenance and run unchanged scientific code."""
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit-workers',type=int,required=True)
    parser.add_argument('command',nargs=argparse.REMAINDER)
    options=parser.parse_args(argv)
    if not 1<=options.fit_workers<=20:parser.error('--fit-workers must be between1 and20')
    command=options.command
    if command and command[0]=='--':command=command[1:]
    args=base.cli.parser().parse_args(command)
    if args.command!='quantify-full' or args.train_only:
        parser.error('only quantify-full fitting is supported')
    out=Path(args.out).resolve()
    out.parent.mkdir(parents=True,exist_ok=True)
    path=Path(str(out)+'.fitworkers.json')
    with Path(str(out)+'.fitworkers.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        # Validates all original scientific and scheduler identities unchanged.
        original_provenance=base.ensure_provenance(args)
        record={'schema':'bulk2cell.fitworkers.v1','output':str(out),
            'extension_sha256':base._sha256(__file__),
            'original_launcher_sha256':base._sha256(base.__file__),
            'original_scheduler_sha256':base._sha256(original_provenance),
            'scientific_identity':base.pipeline._identity(args),
            'effective_fit_workers':options.fit_workers,'fit_pending_limit':options.fit_workers*2,
            'legacy_cli_workers':args.workers,'numerical_threads_per_worker':1,
            'scope':'fit scheduling only; original training, source, checkpoints and inference unchanged'}
        if path.exists():
            if json.loads(path.read_text())!=record:raise ValueError('fit-worker provenance differs; refusing resume')
        else:base.pipeline._atomic_json(path,record)
        original_context=base.patched_scheduler
        base.patched_scheduler=lambda:expanded_scheduler(options.fit_workers)
        try:status=base.main(['--']+command)
        finally:base.patched_scheduler=original_context
        if status==0:
            for directory in [out]+[out/v for v in base.pipeline.VARIANTS]:
                manifest=directory/'run.json'
                data=json.loads(manifest.read_text())
                data['fit_worker_provenance']={'path':str(path),'sha256':base._sha256(path)}
                data['execution_controls'].update(workers=options.fit_workers,
                    fit_pending_limit=options.fit_workers*2,legacy_cli_workers=args.workers)
                base.pipeline._atomic_json(manifest,data)
        return status


if __name__=='__main__':
    raise SystemExit(main())
