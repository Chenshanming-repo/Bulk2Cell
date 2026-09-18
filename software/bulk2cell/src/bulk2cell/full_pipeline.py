"""Bounded, resumable gene fits with a single globally trained capture model.

The training pass follows catalog gene order and sorted molecule order exactly
as the scalar pipeline. SQLite checkpoints hold only the reservoir and sparse
completed gene results. BAM evidence is fetched again in the fit pass, so no
whole-sample read/molecule collection or dense output matrix is materialized.
Linux fork workers share immutable annotation; at most 128 training tasks or
twice the worker count of fitting tasks can be outstanding. A failed gene aborts completion.
"""
from collections import Counter, defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path
import csv
import gzip
import hashlib
import json
import multiprocessing
import os
import pickle
import resource
import sqlite3
import time
import zlib
import numpy as np
from . import adapters, __version__
from .inference import compatible, distance, fit_distance, likelihoods, quantify_gene
from .pipeline import fingerprint, write_tsv

_CONTEXT = {}
_BAM = None
VARIANTS = ('empirical_tes', 'annotation_tes', 'sr_only')


def _initialize_worker():
    """Open one indexed BAM per worker and cap numerical library threading."""
    global _BAM
    import pysam
    from threadpoolctl import threadpool_limits
    threadpool_limits(limits=1)
    _BAM = pysam.AlignmentFile(str(_CONTEXT['args'].bam), 'rb')


def _evidence(gid):
    """Extract one gene while retaining whole-catalog GX ambiguity rules."""
    ctx = _CONTEXT
    molecules, qc = adapters.extract_molecules(ctx['args'].bam, ctx['models'][gid], {gid}, ctx['barcode_set'],
        gene_normalization=ctx['normalization'], bam_handle=_BAM)
    if ctx['args'].max_likelihood_entries > 0 and len(molecules) * len(ctx['models'][gid]) > ctx['args'].max_likelihood_entries:
        raise ValueError(f'{gid}: likelihood budget exceeded ({len(molecules)} molecules, {len(ctx["models"][gid])} models)')
    return molecules, qc


def _train_gene(gid):
    """Return structurally unique distances in scalar molecule iteration order."""
    started = time.monotonic()
    molecules, qc = _evidence(gid)
    ts = _CONTEXT['models'][gid]
    training = []
    cache = {}
    for molecule in molecules:
        key = (molecule.blocks, molecule.junctions)
        if key in cache:
            if cache[key] is not None:
                training.append(cache[key])
            continue
        hit = None
        for transcript in ts:
            if compatible(molecule, transcript):
                if hit is not None:
                    hit = None
                    break
                hit = transcript
        value = distance(molecule, hit) if hit is not None else None
        cache[key] = value
        if value is not None:
            training.append(value)
    return gid, training, qc, dict(worker_seconds=time.monotonic()-started,
        peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)


def _compact_result(result):
    """Replace per-gene dense count arrays with sparse coordinate triples."""
    for key in ('isoform_counts', 'group_counts'):
        matrix = result[key]
        rows, cols = np.nonzero(matrix)
        result[key] = (rows, cols, matrix[rows, cols])
    return result


def cached_likelihoods(molecules, transcripts, model, tes):
    """Evaluate each distinct evidence shape once, preserving scalar arithmetic.

    Cell and UMI identifiers never enter compatibility or distance. Reindexing
    identical rows after the original likelihood kernel is therefore exact,
    including row normalization and empirical TES marginalization.
    """
    representatives = []
    indices = {}
    inverse = []
    for molecule in molecules:
        key = (molecule.gene_id, molecule.blocks, molecule.junctions)
        if key not in indices:
            indices[key] = len(representatives)
            representatives.append(molecule)
        inverse.append(indices[key])
    matrix = likelihoods(representatives, transcripts, model, tes, normalize=True)
    if len(representatives) == len(molecules):
        return matrix
    return matrix[np.asarray(inverse, dtype=np.intp)]


def batched_em(matrices, prior, strength, tolerance, max_iter):
    """Vectorize equal-row cells with exactly the scalar summation axes/order.

    Each cell freezes at its own original convergence iteration. Bucketing by
    row count avoids padding, changed summation order, or cross-cell coupling.
    The caller limits the overall batch; no whole-sample tensor is allocated.
    """
    buckets = defaultdict(list)
    for i, matrix in enumerate(matrices):
        buckets[len(matrix)].append(i)
    answers = [None] * len(matrices)
    for nrows, indices in buckets.items():
        values = np.stack([matrices[i] for i in indices])
        abundance = np.full((len(indices), len(prior)), 1 / len(prior))
        active = np.ones(len(indices), dtype=bool)
        iterations = np.full(len(indices), max_iter, dtype=int)
        converged = np.zeros(len(indices), dtype=bool)
        for iteration in range(1, max_iter + 1):
            where = np.flatnonzero(active)
            posterior = values[where] * abundance[where, None, :]
            posterior /= posterior.sum(axis=2, keepdims=True)
            counts = posterior.sum(axis=1)
            updated = (counts + strength * prior) / (nrows + strength)
            updated = np.maximum(updated, 1e-300)
            updated /= updated.sum(axis=1, keepdims=True)
            done = np.max(np.abs(updated - abundance[where]), axis=1) < tolerance
            abundance[where] = updated
            finished = where[done]
            active[finished] = False
            converged[finished] = True
            iterations[finished] = iteration
            if not active.any():
                break
        posterior = values * abundance[:, None, :]
        posterior /= posterior.sum(axis=2, keepdims=True)
        counts = posterior.sum(axis=1)
        for j, original in enumerate(indices):
            answers[original] = (counts[j], bool(converged[j]), int(iterations[j]))
    return answers


def _fit_gene(gid):
    """Fit three variants while sharing evidence and the annotation likelihood."""
    ctx = _CONTEXT
    args = ctx['args']
    started = time.monotonic()
    molecules, qc = _evidence(gid)
    extraction_seconds = time.monotonic()-started
    ts = sorted(ctx['models'][gid], key=lambda t: t.id)
    common = dict(tau=args.tau, window=args.window, model=ctx['model'],
                  em_strength=args.em_strength, max_iter=args.max_iter, em_batch=batched_em)
    started = time.monotonic()
    annotation = cached_likelihoods(molecules, ts, ctx['model'], {})
    shared_likelihood_seconds = time.monotonic()-started
    started = time.monotonic()
    annotation_fit = _compact_result(quantify_gene(molecules, ts, ctx['sr'], {},
        precomputed_likelihood=annotation, **common))
    annotation_seconds = time.monotonic()-started
    started = time.monotonic()
    sr_fit = _compact_result(quantify_gene(molecules, [replace(t, lr_count=0) for t in ts], ctx['sr'], {},
        precomputed_likelihood=annotation, **common))
    sr_seconds = time.monotonic()-started
    del annotation
    started = time.monotonic()
    if any(t.id in ctx['tes'] for t in ts):
        empirical_matrix = cached_likelihoods(molecules, ts, ctx['model'], ctx['tes'])
        empirical = _compact_result(quantify_gene(molecules, ts, ctx['sr'], ctx['tes'],
            precomputed_likelihood=empirical_matrix, **common))
    else:
        empirical = annotation_fit
    metrics = dict(extraction_seconds=extraction_seconds, shared_annotation_likelihood_seconds=shared_likelihood_seconds,
        variant_worker_seconds=dict(empirical_tes=time.monotonic()-started, annotation_tes=annotation_seconds, sr_only=sr_seconds),
        peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    return gid, dict(empirical_tes=empirical, annotation_tes=annotation_fit, sr_only=sr_fit), len(molecules), qc, metrics


def _ordered_results(function, genes, workers, pending_limit=None):
    """Yield in input order with bounded lookahead and fail on worker errors."""
    pending_limit = workers * 2 if pending_limit is None else pending_limit
    if workers == 1:
        _initialize_worker()
        try:
            for gid in genes:
                yield function(gid)
        finally:
            _BAM.close()
        return
    with ProcessPoolExecutor(max_workers=workers, mp_context=multiprocessing.get_context('fork'),
                             initializer=_initialize_worker) as pool:
        pending = deque()
        iterator = iter(genes)
        for _ in range(pending_limit):
            gid = next(iterator, None)
            if gid is not None:
                pending.append(pool.submit(function, gid))
        while pending:
            yield pending.popleft().result()
            gid = next(iterator, None)
            if gid is not None:
                pending.append(pool.submit(function, gid))


def _pack(value):
    """Compress a trusted local checkpoint payload for SQLite storage."""
    return zlib.compress(pickle.dumps(value, protocol=5), level=1)


def _unpack(value):
    """Read only task-owned local checkpoint data, never external pickle inputs."""
    return pickle.loads(zlib.decompress(value))


def _atomic_json(path, value):
    """Publish a manifest only after fully writing a temporary sibling."""
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, allow_nan=False))
    os.replace(temporary, path)


def _identity(args):
    """Bind resumes to source metadata, small-input hashes, code and parameters."""
    parameters = {k: v for k, v in vars(args).items() if k not in ('out', 'resume', 'workers', 'max_likelihood_entries')}
    paths = [getattr(args, k, None) for k in ('catalog', 'reference', 'isoseq_gff', 'classification', 'bam', 'barcodes', 'salmon', 'tes', 'training_cache')]
    for candidate in (Path(str(args.bam) + '.bai'), Path(args.bam).with_suffix('.bai'), Path(str(args.bam) + '.csi')):
        if candidate.exists():
            paths.append(str(candidate))
    code = {name: hashlib.sha256(Path(__file__).with_name(name).read_bytes()).hexdigest()
            for name in ('full_pipeline.py', 'inference.py', 'adapters.py', 'models.py')}
    return dict(parameters=parameters, inputs=[fingerprint(p) for p in paths if p], code=code)


def _training_identity(args):
    """Bind capture training only to evidence, order, seed, and compatible code."""
    record = _identity(args)
    sources = {str(Path(getattr(args, k)).resolve()) for k in
        ('catalog', 'reference', 'isoseq_gff', 'classification', 'bam', 'barcodes') if getattr(args, k, None)}
    sources.update(str(p.resolve()) for p in
        (Path(str(args.bam)+'.bai'), Path(args.bam).with_suffix('.bai'), Path(str(args.bam)+'.csi')) if p.exists())
    return dict(inputs=[p for p in record['inputs'] if p['path'] in sources],
        parameters={k: getattr(args, k, None) for k in ('genes', 'seed')}, code=record['code'])


def _merge_variant(out, variant, db, models, barcodes, training, manifest):
    """Stream deterministic TSVs and sparse MatrixMarket triples without global lists."""
    dest = out / variant
    dest.mkdir(exist_ok=True)
    headers = {
        'isoforms': ['transcript_id', 'gene_id', 'source', 'lr_count', 'hybrid_prior'],
        'groups': ['group_id', 'gene_id', 'members', 'interpretation'],
        'membership': ['group_id', 'transcript_id', 'within_group_weight', 'interpretation'],
        'gene_qc': ['gene_id', 'isoforms', 'groups', 'input_molecules', 'assigned_molecules', 'incompatible_molecules', 'em_nonconverged_cells', 'max_em_iterations']}
    handles = {k: (dest / (k + '.tsv')).open('w', newline='') for k in headers}
    writers = {k: csv.writer(f, delimiter='\t', lineterminator='\n') for k, f in handles.items()}
    for k in headers:
        writers[k].writerow(headers[k])
    bodies = {k: (dest / (k + '.body')).open('w') for k in ('isoform_counts', 'group_counts')}
    nnz = Counter(); sums = Counter(); totals = Counter(); offsets = dict(isoform_counts=0, group_counts=0)
    cell_index = {c: i for i, c in enumerate(barcodes)}
    try:
        for gid in sorted(models):
            row = db.execute('SELECT result FROM fits WHERE gene=?', (gid,)).fetchone()
            if row is None:
                raise RuntimeError(f'missing completed gene: {gid}')
            _, results, count, _, _ = _unpack(row[0]); result = results[variant]
            ts = sorted(models[gid], key=lambda t: t.id)
            for j, t in enumerate(ts):
                writers['isoforms'].writerow((t.id, gid, t.source, 0 if variant == 'sr_only' else t.lr_count, float(result['prior'][j])))
            for j, indices in enumerate(result['groups']):
                group_id = f'{gid}:group{j+1}'
                label = 'prior_informed_decomposition' if len(indices) > 1 else 'single_member_group'
                writers['groups'].writerow((group_id, gid, len(indices), label))
                denominator = result['prior'][indices].sum()
                for k in indices:
                    weight = float(result['prior'][k] / denominator) if denominator else 1 / len(indices)
                    writers['membership'].writerow((group_id, ts[k].id, weight, label))
            for key in bodies:
                rows, cols, values = result[key]
                for r, c, v in zip(rows, cols, values):
                    bodies[key].write(f'{cell_index[result["cells"][r]]+1} {int(c)+offsets[key]+1} {float(v):.17g}\n')
                nnz[key] += len(values); sums[key] += float(values.sum())
            offsets['isoform_counts'] += len(ts); offsets['group_counts'] += len(result['groups'])
            qc = result['qc']
            writers['gene_qc'].writerow((gid, len(ts), len(result['groups']), count, qc['assigned_molecules'],
                qc['incompatible_molecules'], qc['em_nonconverged_cells'], qc['max_em_iterations']))
            for key in ('assigned_molecules', 'incompatible_molecules', 'em_nonconverged_cells'):
                totals[key] += qc[key]
    finally:
        for handle in list(handles.values()) + list(bodies.values()):
            handle.close()
    for key in bodies:
        if not np.isclose(sums[key], totals['assigned_molecules']):
            raise RuntimeError(f'{variant}: molecule conservation failed')
        target = dest / (key + '.mtx.gz')
        with gzip.open(str(target) + '.tmp', 'wt') as handle:
            handle.write(f'%%MatrixMarket matrix coordinate real general\n% cells x features; fractional expected UMI counts\n{len(barcodes)} {offsets[key]} {nnz[key]}\n')
            with (dest / (key + '.body')).open() as source:
                for block in iter(lambda: source.read(1024 * 1024), ''):
                    handle.write(block)
        os.replace(str(target) + '.tmp', target)
        (dest / (key + '.body')).unlink()
    write_tsv(dest / 'barcodes.tsv', ['barcode'], ((c,) for c in barcodes))
    write_tsv(dest / 'distance_training.tsv', ['distance'], ((d,) for d in training))
    with (dest / 'transcripts.gtf').open('w') as handle:
        for ts in models.values():
            for t in ts:
                for a, b in t.exons:
                    handle.write(f'{t.chrom}\tbulk2cell\texon\t{a+1}\t{b}\t.\t{t.strand}\t.\tgene_id "{t.gene_id}"; transcript_id "{t.id}";\n')
    result = dict(manifest, status='complete', variant=variant, cells=len(barcodes), genes=len(models),
                  isoforms=offsets['isoform_counts'], groups=offsets['group_counts'], qc=dict(totals))
    _atomic_json(dest / 'run.json', result)
    return result


def _run_full_quantification(args):
    """Run or safely resume all three methods using one shared global reservoir."""
    global _CONTEXT
    started = time.monotonic()
    if args.workers < 1:
        raise ValueError('--workers must be positive')
    if args.train_only and (args.salmon or args.tes or args.training_cache):
        raise ValueError('--train-only requires evidence inputs only, without Salmon, TES, or another training cache')
    identity = _training_identity(args) if args.train_only else _identity(args)
    out = Path(args.out)
    if out.exists():
        if not args.resume:
            raise FileExistsError(f'output already exists: {out}; use --resume')
        if not (out / 'identity.json').exists() or json.loads((out / 'identity.json').read_text()) != identity:
            raise ValueError('resume input fingerprint or parameters differ')
    else:
        out.mkdir(parents=True)
        _atomic_json(out / 'identity.json', identity)
    for marker in [out / 'run.json'] + [out / v / 'run.json' for v in VARIANTS]:
        marker.unlink(missing_ok=True)
    if args.catalog:
        transcripts, annotation_qc = adapters.load_regionquant_catalog(args.catalog)
    else:
        transcripts, annotation_qc = adapters.load_transcript_catalog(args.reference, args.isoseq_gff, args.classification)
    gene_aliases = annotation_qc.get('gene_aliases', {})
    selected = {gene_aliases.get(g, g) for g in args.genes} if args.genes else {t.gene_id for t in transcripts}
    if selected - {t.gene_id for t in transcripts}:
        raise ValueError('unknown selected genes')
    transcripts = [t for t in transcripts if t.gene_id in selected]
    if not transcripts:
        raise ValueError('no transcripts selected')
    models = defaultdict(list); normalization = defaultdict(set)
    for t in transcripts:
        models[t.gene_id].append(t)
        normalization[adapters.normalize_gene_id(t.gene_id)].add(t.gene_id)
    barcodes = sorted(adapters.load_barcodes(args.barcodes))
    if not barcodes:
        raise ValueError('barcode whitelist is empty')
    aliases = annotation_qc.get('aliases', {})
    sr, sr_qc = adapters.load_salmon_quant(args.salmon, transcripts, aliases=aliases) if args.salmon else ({}, {'fallback': 'uniform_within_gene'})
    tes, tes_qc = adapters.load_tes_tsv(args.tes) if args.tes else ({}, {'fallback': 'annotation_TES'})
    canonical_tes = defaultdict(list); ids = {t.id for t in transcripts}
    for key, points in tes.items():
        target = aliases.get(key, key)
        if target in ids:
            canonical_tes[target].extend(points)
    tes_qc['matched_transcripts'] = len(canonical_tes)
    _CONTEXT = dict(args=args, models=models, normalization=normalization, barcode_set=set(barcodes), sr=sr, tes=canonical_tes)
    db = sqlite3.connect(out / 'checkpoints.sqlite')
    db.execute('PRAGMA journal_mode=WAL')
    db.execute('PRAGMA synchronous=NORMAL')
    db.execute('CREATE TABLE IF NOT EXISTS state (key TEXT PRIMARY KEY, value BLOB NOT NULL)')
    db.execute('CREATE TABLE IF NOT EXISTS fits (gene TEXT PRIMARY KEY, result BLOB NOT NULL, metrics TEXT NOT NULL)')
    previous = db.execute("SELECT value FROM state WHERE key='elapsed'").fetchone()
    previous_elapsed = _unpack(previous[0]) if previous else 0.0
    def checkpoint_time():
        """Persist observed cumulative wall time at each durable gene checkpoint."""
        elapsed = previous_elapsed + time.monotonic()-started
        db.execute("INSERT OR REPLACE INTO state VALUES ('elapsed',?)", (_pack(elapsed),))
        return elapsed
    try:
        state_row = db.execute("SELECT value FROM state WHERE key='training'").fetchone()
        rng = np.random.default_rng(args.seed)
        completed = 0; seen = 0; training = []; bam_qc = Counter()
        training_metrics = dict(worker_seconds=0.0, peak_rss_mib=0.0)
        if state_row:
            completed, seen, training, rng_state, bam_qc, training_metrics = _unpack(state_row[0]); rng.bit_generator.state = rng_state
        external_training_elapsed = 0.0
        if args.training_cache and not state_row:
            cached = json.loads(Path(args.training_cache).read_text())
            if cached.get('status') != 'training_complete' or cached.get('identity') != _training_identity(args):
                raise ValueError('training cache identity differs from evidence or training parameters')
            completed = len(models); seen = cached['seen']; training = cached['training']
            bam_qc = Counter(cached['bam_qc']); training_metrics = cached['training_metrics']
            rng.bit_generator.state = cached['rng_state']
            external_training_elapsed = cached['elapsed_seconds']
            db.execute("INSERT OR REPLACE INTO state VALUES ('training',?)",
                (_pack((completed, seen, training, rng.bit_generator.state, bam_qc, training_metrics)),))
            db.commit()
        elif args.training_cache:
            external_training_elapsed = json.loads(Path(args.training_cache).read_text())['elapsed_seconds']
        gene_order = list(models)
        for gid, distances, qc, metrics in _ordered_results(_train_gene, gene_order[completed:], args.workers, pending_limit=128):
            for d in distances:
                seen += 1
                if len(training) < 2000:
                    training.append(d)
                else:
                    j = int(rng.integers(seen))
                    if j < 2000:
                        training[j] = d
            completed += 1; bam_qc.update(qc)
            training_metrics['worker_seconds'] += metrics['worker_seconds']
            training_metrics['peak_rss_mib'] = max(training_metrics['peak_rss_mib'], metrics['peak_rss_mib'])
            db.execute("INSERT OR REPLACE INTO state VALUES ('training',?)", (_pack((completed, seen, training, rng.bit_generator.state, bam_qc, training_metrics)),))
            checkpoint_time(); db.commit()
            if completed % 100 == 0:
                _atomic_json(out / 'progress.json', dict(stage='training', completed_genes=completed, genes=len(models)))
        if args.train_only:
            result = dict(status='training_complete', identity=identity, training=training, seen=seen,
                execution_controls=dict(workers=args.workers, training_pending_limit=128),
                rng_state=rng.bit_generator.state, bam_qc=dict(bam_qc), training_metrics=training_metrics,
                elapsed_seconds=checkpoint_time(), cells=len(barcodes), genes=len(models), isoforms=len(transcripts), groups=0, qc={})
            db.commit()
            _atomic_json(out / 'training.json', result)
            return result
        _CONTEXT['model'] = fit_distance(training, args.bandwidth)
        finished = {row[0] for row in db.execute('SELECT gene FROM fits')}
        for result in _ordered_results(_fit_gene, (g for g in sorted(models) if g not in finished), args.workers):
            db.execute('INSERT INTO fits VALUES (?,?,?)', (result[0], _pack(result), json.dumps(result[4])))
            checkpoint_time(); db.commit(); finished.add(result[0])
            if len(finished) % 100 == 0:
                _atomic_json(out / 'progress.json', dict(stage='fitting', completed_genes=len(finished), genes=len(models)))
        if finished != set(models):
            raise RuntimeError('gene checkpoint set is incomplete or unexpected')
        fit_metrics = dict(extraction_seconds=0.0, shared_annotation_likelihood_seconds=0.0,
            variant_worker_seconds={v:0.0 for v in VARIANTS}, peak_rss_mib=0.0)
        for (encoded,) in db.execute('SELECT metrics FROM fits'):
            metrics = json.loads(encoded)
            for key in ('extraction_seconds', 'shared_annotation_likelihood_seconds'):
                fit_metrics[key] += metrics[key]
            for variant in VARIANTS:
                fit_metrics['variant_worker_seconds'][variant] += metrics['variant_worker_seconds'][variant]
            fit_metrics['peak_rss_mib'] = max(fit_metrics['peak_rss_mib'], metrics['peak_rss_mib'])
        manifest = dict(software='bulk2cell', version=__version__, matrix_orientation='cells_by_features',
            parameters=vars(args), execution_controls=dict(workers=args.workers, training_pending_limit=128, fit_pending_limit=args.workers*2),
            inputs=identity['inputs'], annotation_qc=annotation_qc, salmon_qc=sr_qc,
            tes_qc=tes_qc, bam_qc=dict(bam_qc), distance_model=dict(kind='reflected_gaussian_kde' if training else 'exponential_fallback',
            bandwidth=args.bandwidth, training_eligible_molecules=seen, training_retained=len(training)),
            runtime_scope='shared_three_variant_runner', external_training_elapsed_seconds=external_training_elapsed, training_metrics=training_metrics, fit_metrics=fit_metrics,
            elapsed_seconds=checkpoint_time(), timing_policy='Cumulative observed wall time through checkpoints; interrupted work after the last checkpoint may be unrecorded', peak_rss_scope='parent_only; worker lifetime maximum reported separately, not aggregate process-tree RSS', peak_rss_mib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
            limitations=['GX single-gene assignments required', 'Within-group decomposition uses fixed bulk priors',
                         'Peak RSS reports parent process only; worker memory is additional'])
        variants = {v: _merge_variant(out, v, db, models, barcodes, training, manifest) for v in VARIANTS}
        elapsed = checkpoint_time(); db.commit()
        for v, record in variants.items():
            record['elapsed_seconds'] = elapsed
            record['peak_rss_mib'] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024
            record['variant_configuration'] = dict(lr_counts='zero' if v == 'sr_only' else 'catalog',
                tes='empirical_with_annotation_fallback' if v == 'empirical_tes' else 'annotation')
            _atomic_json(out / v / 'run.json', record)
        result = dict(variants['empirical_tes'], variants={v: str(out / v) for v in VARIANTS},
                      elapsed_seconds=elapsed)
        _atomic_json(out / 'run.json', result)
        _atomic_json(out / 'progress.json', dict(stage='complete', completed_genes=len(models), genes=len(models)))
        return result
    finally:
        db.close()


def run_full_quantification(args):
    """Hold an exclusive advisory lock throughout checkpoint and output mutation."""
    import fcntl
    out = Path(args.out).resolve()
    out.parent.mkdir(parents=True, exist_ok=True)
    with (out.parent / ('.' + out.name + '.lock')).open('a') as lock:
        try:
            fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError as error:
            raise RuntimeError(f'quantification already running for {out}') from error
        return _run_full_quantification(args)
