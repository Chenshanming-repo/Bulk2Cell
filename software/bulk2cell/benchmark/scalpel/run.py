#!/usr/bin/env python3
"""Run pinned upstream SCALPEL with resumable chromosome-bounded stages."""
from __future__ import annotations

import argparse
from bisect import bisect_right
from collections import Counter, OrderedDict, defaultdict
from concurrent.futures import as_completed, ThreadPoolExecutor
import csv
import fcntl
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import shutil
import subprocess
import signal
import threading
import time

import pysam

PIN = "3c31ffc6aa3b7e422623aececfaea32b3aa358ca"
PARAMETERS = {"distance_transcriptome": 600, "distance_exon": 30,
              "gene_fraction": "98%", "binsize": 20, "ip_distance": 60}


def bounded_map(function, items, workers):
    """Run independent tasks with a fail-fast barrier and bounded concurrency."""
    started = time.monotonic()
    items = list(items)
    if workers == 1:
        return [function(item) for item in items], time.monotonic() - started
    cancelled = threading.Event()
    with _ACTIVE_CHILDREN_LOCK:
        if _ACTIVE_CHILDREN:
            raise RuntimeError("Cannot start a stage barrier with active child processes")
        _GLOBAL_CANCEL.clear()

    def invoke(item):
        """Stop the stage at its worker boundary after the first task failure."""
        if cancelled.is_set():
            raise RuntimeError("Stage cancelled after another worker failed")
        try:
            return function(item)
        except BaseException:
            cancelled.set()
            _GLOBAL_CANCEL.set()
            terminate_active_children()
            raise

    executor = ThreadPoolExecutor(max_workers=workers)
    futures = {}
    results = [None] * len(items)
    try:
        for index, item in enumerate(items):
            futures[executor.submit(invoke, item)] = index
        for future in as_completed(futures):
            results[futures[future]] = future.result()
    except BaseException:
        cancelled.set()
        _GLOBAL_CANCEL.set()
        for future in futures:
            future.cancel()
        terminate_active_children()
        executor.shutdown(wait=True, cancel_futures=True)
        raise
    executor.shutdown(wait=True)
    return results, time.monotonic() - started


_ACTIVE_CHILDREN = set()
_ACTIVE_CHILDREN_LOCK = threading.Lock()
_GLOBAL_CANCEL = threading.Event()
_MAIN_LAUNCH_CRITICAL = False
_SIGNAL_CANCELLED = False


def terminate_active_children(wait=True):
    """Terminate owned process groups and optionally wait before forced killing."""
    with _ACTIVE_CHILDREN_LOCK:
        processes = list(_ACTIVE_CHILDREN)
    for process in processes:
        if process.poll() is None:
            os.killpg(process.pid, signal.SIGTERM)
    if not wait:
        return
    deadline = time.monotonic() + 5
    for process in processes:
        remaining = max(0, deadline - time.monotonic())
        try:
            process.wait(timeout=remaining)
        except subprocess.TimeoutExpired:
            if process.poll() is None:
                os.killpg(process.pid, signal.SIGKILL)


def install_signal_cleanup():
    """Install parent signal handlers that terminate and reap owned R groups."""
    def handle(signum, frame):
        """Set a lock-free flag and unwind after atomic main-thread registration."""
        global _SIGNAL_CANCELLED
        _SIGNAL_CANCELLED = True
        if _MAIN_LAUNCH_CRITICAL:
            return
        raise SystemExit(128 + signum)
    signal.signal(signal.SIGTERM, handle)
    signal.signal(signal.SIGINT, handle)


def exon_rows(path):
    """Yield one-based closed GTF exons with stable gene and transcript IDs."""
    with Path(path).open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "exon":
                continue
            attrs = dict(re.findall(r'(\w+) "([^"]*)"', fields[8]))
            yield (fields[0], int(fields[3]), int(fields[4]), fields[6],
                   attrs["gene_id"], attrs["transcript_id"])


def file_fingerprint(path):
    """Return a content fingerprint suitable for binding resumable stages."""
    path = Path(path).resolve()
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    stat = path.stat()
    return {"path": str(path), "size": stat.st_size, "sha256": digest.hexdigest()}


def validate_prepared_reads(bundle, bam, barcodes, required_contigs):
    """Validate an atomic shared BED certificate and return required BED paths."""
    bundle = Path(bundle).resolve()
    certificate_path = bundle/"certificate.json"
    if not certificate_path.is_file():
        raise FileNotFoundError(f"Prepared-read certificate is missing: {certificate_path}")
    certificate = json.loads(certificate_path.read_text())
    if certificate.get("schema") != "bulk2cell.shared_reads.v1" or not certificate.get("completed"):
        raise ValueError("Prepared-read certificate is incomplete or has an unknown schema")
    if certificate.get("cross_cb_read_ids") != 0:
        raise ValueError("Prepared reads contain read IDs under multiple cell barcodes")
    expected_rules = {"mapped", "primary", "not supplementary", "not QC-fail",
                      "CB in exact 10,250-cell whitelist", "nonempty UB"}
    if set(certificate.get("selection_rules", [])) != expected_rules:
        raise ValueError("Prepared-read selection rules do not match the harness")
    if certificate.get("read_id_rule") != "query_name before first /":
        raise ValueError("Prepared-read ID normalization does not match upstream")
    if int(certificate.get("projection_records", 0)) < 1:
        raise ValueError("Prepared-read certificate has no projected alignments")
    def matches(path, record):
        """Check one producer record against current file content."""
        observed = file_fingerprint(path)
        return (observed["size"] == record.get("bytes") and
                observed["sha256"] == record.get("sha256"))
    if not matches(__file__, certificate["scalpel_consumer_source"]):
        raise ValueError("Prepared-read consumer source fingerprint does not match")
    if not matches(bam, certificate["bam"]):
        raise ValueError("Prepared-read BAM fingerprint does not match")
    if not matches(barcodes, certificate["canonical_sorted_barcodes"]):
        raise ValueError("Prepared-read barcode fingerprint does not match")
    bed_paths = {}
    for contig in sorted(required_contigs):
        filename = contig.replace("/", "_") + ".bed"
        if filename not in certificate.get("bed_files", {}):
            raise ValueError(f"Prepared reads are missing contig: {contig}")
        path = bundle/"bed"/filename
        if not matches(path, certificate["bed_files"][filename]):
            raise ValueError(f"Prepared BED fingerprint does not match: {contig}")
        bed_paths[contig] = path
    return {"certificate": certificate, "certificate_path": certificate_path,
            "bed_paths": bed_paths}


def stage_identity(stage, sources, parameters):
    """Build a deterministic stage identity from source contents and parameters."""
    return {"stage": stage,
            "sources": {key: file_fingerprint(value) for key, value in sorted(sources.items())},
            "parameters": parameters}


def write_stage_manifest(path, identity, outputs, metrics):
    """Atomically publish a completed stage and its verified output inventory."""
    path = Path(path)
    resolved_outputs = [Path(item).resolve() for item in outputs]
    payload = {**identity, "complete": True,
               "outputs": [str(item) for item in resolved_outputs],
               "output_fingerprints": [file_fingerprint(item)
                                       for item in resolved_outputs],
               "metrics": metrics}
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True))
    temporary.replace(path)


def completed_stage(path, identity):
    """Return a valid completed manifest or reject stale/corrupt resume state."""
    path = Path(path)
    if not path.exists():
        return None
    payload = json.loads(path.read_text())
    observed = {key: payload.get(key) for key in ("stage", "sources", "parameters")}
    if observed != identity:
        raise ValueError(f"Resume manifest {path} does not match current sources/parameters")
    if not payload.get("complete"):
        raise ValueError(f"Resume manifest is not complete: {path}")
    outputs = payload.get("outputs", [])
    fingerprints = payload.get("output_fingerprints")
    if not isinstance(fingerprints, list) or len(fingerprints) != len(outputs):
        raise ValueError(f"Resume manifest lacks exact output fingerprints: {path}")
    for output, expected in zip(outputs, fingerprints):
        output = Path(output)
        if not output.is_file():
            raise FileNotFoundError(f"Resume manifest declared output is missing: {output}")
        if file_fingerprint(output) != expected:
            raise ValueError(f"Resume output fingerprint changed: {output}")
    return payload


def acquire_run_lock(directory):
    """Acquire a nonblocking exclusive writer lock for one output directory."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    handle = (directory/".run.lock").open("a")
    try:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as error:
        handle.close()
        raise RuntimeError(f"Another SCALPEL writer is already active: {directory}") from error
    return handle


class HandlePool:
    """Keep a bounded number of append-only text shard handles open."""

    def __init__(self, limit=32):
        """Initialize an LRU pool with at most `limit` open files."""
        self.limit = limit
        self.handles = OrderedDict()

    def writer(self, path):
        """Return a TSV writer while enforcing the open-file limit."""
        path = Path(path)
        if path in self.handles:
            handle = self.handles.pop(path)
        else:
            if len(self.handles) >= self.limit:
                _, old = self.handles.popitem(last=False)
                old.close()
            handle = path.open("a", newline="")
        self.handles[path] = handle
        return csv.writer(handle, delimiter="\t", lineterminator="\n")

    def close(self):
        """Close every pooled handle."""
        for handle in self.handles.values():
            handle.close()
        self.handles.clear()


def assign_cb_batches(record_counts, max_cells_per_batch):
    """Assign each cell once using deterministic load-balanced bounded batches."""
    if max_cells_per_batch < 1:
        raise ValueError("max_cells_per_batch must be positive")
    cells = sorted(record_counts, key=lambda cell: (-record_counts[cell], cell))
    batch_count = max(1, (len(cells) + max_cells_per_batch - 1) // max_cells_per_batch)
    loads = [0] * batch_count
    sizes = [0] * batch_count
    result = {}
    for cell in cells:
        choices = [i for i in range(batch_count) if sizes[i] < max_cells_per_batch]
        batch = min(choices, key=lambda i: (loads[i], sizes[i], i))
        result[cell] = batch
        loads[batch] += record_counts[cell]
        sizes[batch] += 1
    return result


def _external_sort_unique(inputs, output, temporary_directory, keys):
    """Merge tabular inputs through GNU sort with bounded process memory."""
    output = Path(output)
    temporary = output.with_suffix(output.suffix + ".tmp")
    command = ["sort", "-u", "-T", str(temporary_directory), "-S", "1G", *keys,
               *[str(Path(path)) for path in inputs], "-o", str(temporary)]
    subprocess.run(command, check=True)
    temporary.replace(output)


def validate_readid_cb_pairs(path, temporary_directory):
    """Fail if any upstream readID is observed under more than one cell barcode."""
    sorted_path = Path(path).with_suffix(".sorted.tsv")
    _external_sort_unique([path], sorted_path, temporary_directory, ["-k1,1", "-k2,2"])
    previous_id = previous_cb = None
    pairs = 0
    with sorted_path.open() as handle:
        for line in handle:
            read_id, cb = line.rstrip("\n").split("\t")
            pairs += 1
            if read_id == previous_id and cb != previous_cb:
                raise ValueError(f"readID occurs under multiple cell barcodes: {read_id}")
            previous_id, previous_cb = read_id, cb
    return {"readid_cb_pairs": pairs, "cross_cb_readids": 0}


def reconcile_observed_pairs(pair_files, unique_genes_path, temporary_directory):
    """Derive globally single-model genes from all observed batch pairs."""
    pairs_path = Path(unique_genes_path).with_name("observed_pairs.tsv")
    _external_sort_unique(pair_files, pairs_path, temporary_directory, ["-k1,1", "-k2,2"])
    temporary = Path(unique_genes_path).with_suffix(".tmp")
    observed = unique_count = 0
    current_gene, transcripts = None, 0
    with pairs_path.open() as source, temporary.open("w") as target:
        for line in source:
            gene, _ = line.rstrip("\n").split("\t")
            observed += 1
            if current_gene is not None and gene != current_gene:
                if transcripts == 1:
                    target.write(current_gene + "\n"); unique_count += 1
                transcripts = 0
            current_gene = gene; transcripts += 1
        if current_gene is not None and transcripts == 1:
            target.write(current_gene + "\n"); unique_count += 1
    temporary.replace(unique_genes_path)
    return {"observed_pairs": observed, "unique_genes": unique_count}


def write_mapping_bridges(source_path, directory):
    """Generate SHA-bound mapping scripts that expose and preserve an upstream bug."""
    source = Path(source_path).read_text()
    expression = "unspliceds = dplyr::filter(unspliceds, !(ftrs %in% trs.todel))"
    prefix, separator, suffix = source.rpartition(expression)
    if not separator or "#---- concordance unspliced/spliced fragments" not in prefix:
        raise ValueError("Pinned final concordance expression was not found")
    hook = ("fwrite(trs.todel, paste0(args$output_path, "
            "'.concordance_candidates.tsv'), sep='\t', col.names=FALSE)\n")
    instrumented_text = prefix + hook + expression + suffix
    noop_text = prefix + hook + "unspliceds = unspliceds" + suffix
    directory = Path(directory)
    instrumented = directory/"mapping_filtering_instrumented.R"
    noop = directory/"mapping_filtering_global_noop.R"
    for path, text_value in ((instrumented, instrumented_text), (noop, noop_text)):
        temporary = path.with_suffix(".R.tmp")
        temporary.write_text(text_value); temporary.replace(path)
    return instrumented, noop


def concordance_correction_plan(candidate_paths, global_path, temporary_directory):
    """Plan reruns that reproduce the upstream `%in% data.table` cardinality bug."""
    _external_sort_unique(candidate_paths, global_path, temporary_directory, [])
    with Path(global_path).open() as handle:
        global_count = sum(1 for _ in handle)
    local_counts = []
    for path in candidate_paths:
        with Path(path).open() as handle:
            local_counts.append(sum(1 for line in handle if line.rstrip("\n")))
    reruns = ([index for index, count in enumerate(local_counts) if count == 1]
              if global_count > 1 else [])
    return {"global_candidates": global_count, "local_counts": local_counts,
            "rerun_indexes": reruns}


def filtered_batch_has_evidence(observed_pairs_path):
    """Return whether successful original filtering retained any observed model."""
    return Path(observed_pairs_path).stat().st_size > 0


def install_diagnostic_sink(directory):
    """Discard the unused upstream distance diagnostic with explicit provenance."""
    directory = Path(directory)
    diagnostic = directory/"read_distance_distribution_on_transcriptomic_scope.txt"
    if diagnostic.exists() or diagnostic.is_symlink():
        if not diagnostic.is_symlink() or diagnostic.resolve() != Path("/dev/null"):
            raise FileExistsError(f"Refusing to replace diagnostic path: {diagnostic}")
        diagnostic.unlink()
    diagnostic.symlink_to("/dev/null")
    provenance = directory/"diagnostic_provenance.json"
    temporary = provenance.with_suffix(".json.tmp")
    temporary.write_text(json.dumps({"path": diagnostic.name, "target": "/dev/null",
        "policy": "discard unused diagnostic only; mapping computation unchanged"}, indent=2))
    temporary.replace(provenance)


def _cell_from_fragment(fragment):
    """Extract a cell barcode from the original CB::UB fragment identifier."""
    if not fragment.startswith("CB:Z:") or "::UB:Z:" not in fragment:
        raise ValueError(f"Malformed CB::UB fragment identifier: {fragment}")
    return fragment[5:].split("::UB:Z:", 1)[0]


def upstream_read_id(read_id):
    """Return the readID produced by upstream tidyr separation at first slash."""
    return read_id.split("/", 1)[0]


def repartition_by_cb(out, chromosome_shards, chromosomes, max_cells_per_batch,
                      cross_cb_certificate=None):
    """Repartition chromosome BEDs by complete cells using bounded file handles."""
    out = Path(out)
    pair_path = out/"readid_cb.tsv"
    batch_chromosomes, batches, qc = {}, [], Counter()
    pair_handle = pair_path.open("w") if cross_cb_certificate is None else None
    try:
        for chromosome_shard in chromosome_shards:
            chromosome_directory = out/chromosome_shard
            counts = Counter()
            with (chromosome_directory/"reads.bed").open() as source:
                for line in source:
                    fields = line.rstrip("\n").split("\t")
                    if len(fields) != 6:
                        raise ValueError(f"Malformed BED row in {chromosome_directory}")
                    cell = _cell_from_fragment(fields[5])
                    counts[cell] += 1
                    if pair_handle is not None:
                        pair_handle.write(f"{upstream_read_id(fields[4])}\t{cell}\n")
            assignments = assign_cb_batches(counts, max_cells_per_batch)
            names = {index: f"{chromosome_shard}_cb_{index:04d}"
                     for index in sorted(set(assignments.values()))}
            for index, name in names.items():
                directory = out/name; directory.mkdir(exist_ok=True)
                annotation = directory/"annotation.tsv"
                annotation.symlink_to((chromosome_directory/"annotation.tsv").resolve())
                (directory/"reads.bed").write_text("")
                batches.append(name); batch_chromosomes[name] = chromosome_shard
            pool = HandlePool()
            try:
                with (chromosome_directory/"reads.bed").open() as source:
                    for line in source:
                        fields = line.rstrip("\n").split("\t")
                        cell = _cell_from_fragment(fields[5])
                        pool.writer(out/names[assignments[cell]]/"reads.bed").writerow(fields)
            finally:
                pool.close()
            (chromosome_directory/"reads.bed").unlink()
            qc["cell_batches"] += len(names)
            qc["cells_partitioned"] += len(assignments)
    finally:
        if pair_handle is not None:
            pair_handle.close()
    if cross_cb_certificate is None:
        qc.update(validate_readid_cb_pairs(pair_path, out))
    else:
        qc.update({"readid_cb_pairs": cross_cb_certificate["unique_read_id_cb_pairs"],
                   "cross_cb_readids": cross_cb_certificate["cross_cb_read_ids"],
                   "prepared_cross_cb_certificate": 1})
    (out/"batch_chromosomes.json").write_text(json.dumps({
        "batch_to_chromosome_shard": batch_chromosomes,
        "chromosomes": chromosomes}, indent=2, sort_keys=True))
    return batches, qc


def _blocks(read):
    """Return intron-split zero-based half-open blocks, retaining deletions."""
    blocks = []
    pos = read.reference_start
    start = pos
    for op, length in read.cigartuples or ():
        if op == 3:
            blocks.append((start, pos))
            pos += length
            start = pos
        elif op in (0, 2, 7, 8):
            pos += length
    blocks.append((start, pos))
    return blocks


def prepare(gtf, quant, bam, barcodes, out, prepared_reads=None):
    """Stream annotations and reads into the original chromosome boundaries.

    Every selected alignment on an annotated chromosome is emitted in full.
    SCALPEL's CB::UB fragment families, cross-gene collisions, range filtering,
    and spliced/unspliced checks therefore see exactly the original context.
    """
    out = Path(out)
    gene_ids, transcript_ids, tx_gene, gene_chrom = set(), set(), {}, {}
    chroms = set()
    for chrom, start, end, strand, gene, tx in exon_rows(gtf):
        chroms.add(chrom); gene_ids.add(gene); transcript_ids.add(tx)
        if tx in tx_gene and tx_gene[tx] != gene:
            raise ValueError(f"Transcript occurs in multiple genes: {tx}")
        tx_gene[tx] = gene
        if gene in gene_chrom and gene_chrom[gene] != chrom:
            raise ValueError(f"Gene occurs on multiple chromosomes: {gene}")
        gene_chrom[gene] = chrom
    genes = {value: f"G{i:06d}" for i, value in enumerate(sorted(gene_ids))}
    transcripts = {value: f"T{i:06d}" for i, value in enumerate(sorted(transcript_ids))}
    aliases = {**{alias: stable for stable, alias in genes.items()},
               **{alias: stable for stable, alias in transcripts.items()}}
    (out/"aliases.json").write_text(json.dumps(aliases, indent=2, sort_keys=True))
    chrom_alias = {chrom: f"chrom_{index:03d}" for index, chrom in enumerate(sorted(chroms))}
    for directory_name in chrom_alias.values():
        directory = out/directory_name
        directory.mkdir(exist_ok=True)
        (directory/"annotation.tsv").write_text("")
        (directory/"reads.bed").write_text("")
    pool = HandlePool()
    try:
        for chrom, start, end, strand, gene, tx in exon_rows(gtf):
            pool.writer(out/chrom_alias[chrom]/"annotation.tsv").writerow(
                [chrom, start, end, end-start+1, strand, genes[gene], genes[gene],
                 transcripts[tx], transcripts[tx]])
    finally:
        pool.close()
    with Path(quant).open() as handle:
        tpm = {row["Name"]: float(row["TPM"]) + 1
               for row in csv.DictReader(handle, delimiter="\t")}
    missing = transcript_ids - set(tpm)
    if missing:
        raise ValueError(f"Salmon quantification is missing {len(missing)} annotated transcripts")
    totals = Counter()
    for tx in sorted(transcript_ids):
        totals[tx_gene[tx]] += tpm[tx]
    with (out/"bulk.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["gene_name", "transcript_name", "bulk_TPMperc"])
        for tx in sorted(transcript_ids):
            writer.writerow([genes[tx_gene[tx]], transcripts[tx],
                             tpm[tx]/totals[tx_gene[tx]]])
    if prepared_reads is not None:
        prepared = validate_prepared_reads(prepared_reads, bam, barcodes, chroms)
        for chrom, path in prepared["bed_paths"].items():
            destination = out/chrom_alias[chrom]/"reads.bed"
            destination.unlink()
            destination.symlink_to(path)
        certificate = prepared["certificate"]
        annotated_blocks = sum(row["blocks"] for row in certificate.get("contig_cb_blocks", [])
                               if row["contig"] in chroms)
        qc = Counter({"input_alignments": certificate.get("projection_records", 0),
                      "bed_blocks": annotated_blocks,
                      "prepared_read_bundle": 1,
                      "cross_cb_readids": certificate["cross_cb_read_ids"]})
        ordered = [chrom_alias[chrom] for chrom in sorted(chroms)]
        payload = {"boundary": "original chromosome", "shards": ordered,
                   "chromosomes": {chrom_alias[chrom]: chrom for chrom in sorted(chroms)},
                   "genes": len(genes), "transcripts": len(transcripts),
                   "prepared_reads": str(prepared["certificate_path"])}
        (out/"shards.json").write_text(json.dumps(payload, indent=2))
        return ordered, qc
    allowed = set(Path(barcodes).read_text().splitlines())
    qc, pool = Counter(), HandlePool()
    try:
        with pysam.AlignmentFile(bam, "rb") as source:
            for read in source.fetch(until_eof=True):
                qc["input_alignments"] += 1
                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    qc["nonprimary_or_unmapped"] += 1
                    continue
                if (not read.has_tag("CB") or not read.has_tag("UB") or
                        read.get_tag("CB") not in allowed):
                    qc["missing_or_unselected_tags"] += 1
                    continue
                if read.reference_name not in chrom_alias:
                    qc["outside_annotation"] += 1
                    continue
                blocks = _blocks(read)
                strand = "-" if read.is_reverse else "+"
                fragment = f"CB:Z:{read.get_tag('CB')}::UB:Z:{read.get_tag('UB')}"
                writer = pool.writer(out/chrom_alias[read.reference_name]/"reads.bed")
                for index, (start, end) in enumerate(blocks, 1):
                    writer.writerow([read.reference_name, start, end, strand,
                                     f"{read.query_name}/{index}", fragment])
                qc["selected_alignments"] += 1
                qc["bed_blocks"] += len(blocks)
    finally:
        pool.close()
    ordered = [chrom_alias[chrom] for chrom in sorted(chroms)]
    payload = {"boundary": "original chromosome", "shards": ordered,
               "chromosomes": {chrom_alias[chrom]: chrom for chrom in sorted(chroms)},
               "genes": len(genes), "transcripts": len(transcripts)}
    (out/"shards.json").write_text(json.dumps(payload, indent=2))
    return ordered, qc


def verify_upstream(upstream):
    """Verify the pinned upstream revision and exact scientific source hashes."""
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"],
                                       cwd=upstream, text=True).strip()
    if revision != PIN:
        raise ValueError(f"Expected upstream {PIN}, found {revision}")
    manifest = json.loads(Path(__file__).with_name("upstream-sha256.json").read_text())
    for relative, expected in manifest.items():
        if hashlib.sha256((upstream/relative).read_bytes()).hexdigest() != expected:
            raise ValueError(f"Modified upstream source: {relative}")
    return revision


def run_command(command, cwd, log_path):
    """Run a command and return wall time and incremental child peak RSS."""
    before = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    started = time.monotonic()
    with Path(log_path).open("a") as log:
        status = None
        global _MAIN_LAUNCH_CRITICAL
        try:
            is_main = threading.current_thread() is threading.main_thread()
            if is_main:
                _MAIN_LAUNCH_CRITICAL = True
            with _ACTIVE_CHILDREN_LOCK:
                if _GLOBAL_CANCEL.is_set() or _SIGNAL_CANCELLED:
                    raise RuntimeError("Stage cancelled before subprocess launch")
                status = subprocess.Popen([str(x) for x in command], cwd=cwd,
                                          stdout=log, stderr=log,
                                          start_new_session=True)
                _ACTIVE_CHILDREN.add(status)
            if is_main:
                _MAIN_LAUNCH_CRITICAL = False
                if _GLOBAL_CANCEL.is_set() or _SIGNAL_CANCELLED:
                    raise SystemExit(128 + signal.SIGTERM)
            returncode = status.wait()
        except BaseException:
            _MAIN_LAUNCH_CRITICAL = False
            _GLOBAL_CANCEL.set()
            if status is not None and status.poll() is None:
                os.killpg(status.pid, signal.SIGTERM)
                try:
                    status.wait(timeout=5)
                except subprocess.TimeoutExpired:
                    if status.poll() is None:
                        os.killpg(status.pid, signal.SIGKILL)
                    status.wait()
            raise
        finally:
            _MAIN_LAUNCH_CRITICAL = False
            if status is not None:
                with _ACTIVE_CHILDREN_LOCK:
                    _ACTIVE_CHILDREN.discard(status)
    elapsed = time.monotonic() - started
    after = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    if returncode:
        raise subprocess.CalledProcessError(returncode, command)
    return {"command": [str(x) for x in command], "cwd": str(cwd),
            "seconds": elapsed, "peak_rss_kib": after,
            "peak_rss_increment_kib": max(0, after-before),
            "returncode": returncode}


def concatenate(paths, destination):
    """Atomically concatenate shard files without materializing them in memory."""
    destination = Path(destination)
    temporary = destination.with_suffix(destination.suffix+".tmp")
    with temporary.open("wb") as target:
        for path in paths:
            with Path(path).open("rb") as source:
                shutil.copyfileobj(source, target, length=8*1024*1024)
    temporary.replace(destination)



def start_invocation_ledger(directory, arguments):
    """Atomically append a running invocation and expose prior timing gaps."""
    path = Path(directory)/"invocations.json"
    payload = json.loads(path.read_text()) if path.is_file() else {"invocations": []}
    incomplete = sum(row.get("status") == "running" for row in payload["invocations"])
    identifier = f"{time.time_ns()}-{os.getpid()}"
    payload["invocations"].append({"id": identifier, "status": "running",
        "started_unix": time.time(), "pid": os.getpid(), "arguments": arguments})
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(payload, indent=2)+"\n")
    temporary.replace(path)
    completed_seconds = sum(float(row.get("wall_seconds", 0))
                            for row in payload["invocations"]
                            if row.get("status") == "complete")
    return identifier, completed_seconds, incomplete


def finish_invocation_ledger(directory, identifier, wall_seconds):
    """Atomically mark the current invocation complete and retain all attempts."""
    path = Path(directory)/"invocations.json"
    payload = json.loads(path.read_text())
    matches = [row for row in payload["invocations"] if row["id"] == identifier]
    if len(matches) != 1 or matches[0]["status"] != "running":
        raise RuntimeError("Invocation ledger entry is missing or not running")
    matches[0].update({"status": "complete", "finished_unix": time.time(),
                       "wall_seconds": wall_seconds})
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(payload, indent=2)+"\n")
    temporary.replace(path)


def validate_prepare_adoption(directory):
    """Require the final marker for an explicitly adopted preparation manifest."""
    directory = Path(directory)
    adoption = directory/"prepare_adoption.json"
    if not adoption.exists():
        return None
    marker_path = directory/"PREPARE_ADOPTION_COMPLETE.json"
    if not marker_path.is_file():
        raise RuntimeError("Prepared-stage adoption transaction is incomplete")
    marker = json.loads(marker_path.read_text())
    expected = {
        "adoption": file_fingerprint(adoption),
        "prepare_manifest": file_fingerprint(directory/"prepare.complete.json")}
    if not marker.get("completed") or any(marker.get(key) != value
                                          for key, value in expected.items()):
        raise ValueError("Prepared-stage adoption completion marker is stale")
    return json.loads(adoption.read_text())


def collect_stage_records(directory):
    """Collect completed-stage metrics for resume-stable timing and commands."""
    records = []
    for path in sorted(Path(directory).rglob("*.complete.json")):
        payload = json.loads(path.read_text())
        if payload.get("complete"):
            records.append({"manifest": str(path), "stage": payload["stage"],
                            "metrics": payload.get("metrics", {})})
    return records


def summarize_resources(records):
    """Sum stage wall times and retain the maximum recorded child RSS."""
    seconds, peak, commands = 0.0, 0, []
    for record in records:
        metrics = record["metrics"]
        items = metrics.get("commands", [metrics])
        for item in items:
            seconds += float(item.get("seconds", 0))
            peak = max(peak, int(item.get("peak_rss_kib", 0)))
            if "command" in item:
                commands.append(item)
    return {"stage_wall_seconds": seconds, "peak_child_rss_kib": peak,
            "commands": commands}

def main():
    """Execute pinned stages with global calibration and chromosome-bounded EM."""
    run_started = time.monotonic()
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("gtf", "quant", "bam", "barcodes", "out", "upstream", "rscript"):
        parser.add_argument("--"+name, required=True, type=Path)
    parser.add_argument("--calibration-probabilities", type=Path)
    parser.add_argument("--prepared-reads", type=Path)
    parser.add_argument("--boundary", choices=("chromosome", "cell-batch"),
                        default="chromosome")
    parser.add_argument("--cells-per-batch", type=int, default=128)
    parser.add_argument("--r-workers", type=int, default=1)
    parser.add_argument("--stop-after-prepare", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    if args.r_workers < 1:
        parser.error("--r-workers must be at least 1")
    for key, value in vars(args).items():
        if isinstance(value, Path):
            setattr(args, key, value.resolve())
    revision = verify_upstream(args.upstream)
    if args.out.exists() and not args.resume:
        raise FileExistsError(args.out)
    args.out.mkdir(parents=True, exist_ok=True)
    lock_handle = acquire_run_lock(args.out)
    install_signal_cleanup()
    global _SIGNAL_CANCELLED
    _SIGNAL_CANCELLED = False
    _GLOBAL_CANCEL.clear()
    if args.resume and any(args.out.glob("*.mtx.gz")):
        raise FileExistsError("Cannot resume a run with published matrices")
    invocation_id, prior_invocation_seconds, incomplete_invocations = (
        start_invocation_ledger(args.out, {
            "resume": args.resume, "r_workers": args.r_workers,
            "boundary": args.boundary, "stop_after_prepare": args.stop_after_prepare}))
    if args.resume:
        (args.out/"summary.json").unlink(missing_ok=True)
        (args.out/"matrix_summary.json").unlink(missing_ok=True)
    run_binding = {
        "runner": file_fingerprint(__file__),
        "em_batch": file_fingerprint(Path(__file__).with_name("em_batch.R")),
        "batch_pairs": file_fingerprint(Path(__file__).with_name("batch_pairs.R")),
        "batch_unique": file_fingerprint(Path(__file__).with_name("batch_unique.R")),
        "upstream_manifest": file_fingerprint(
            Path(__file__).with_name("upstream-sha256.json")),
        "upstream_commit": revision}
    def bound(parameters):
        """Bind every resumable stage to harness and upstream source identity."""
        return {**parameters, "run_binding": run_binding}
    inputs = {"gtf": args.gtf, "quant": args.quant, "bam": args.bam,
              "barcodes": args.barcodes}
    prepare_arguments = dict(inputs)
    prepare_arguments["prepared_reads"] = args.prepared_reads
    if args.prepared_reads is not None:
        inputs["prepared_certificate"] = args.prepared_reads/"certificate.json"
    prep_identity = stage_identity("prepare", inputs, bound({
        "boundary": args.boundary, "cells_per_batch": args.cells_per_batch
        if args.boundary == "cell-batch" else None}))
    prep_manifest = args.out/"prepare.complete.json"
    resumed = completed_stage(prep_manifest, prep_identity) if args.resume else None
    if resumed:
        validate_prepare_adoption(args.out)
        shards = json.loads((args.out/"shards.json").read_text())["shards"]
        qc = Counter(resumed["metrics"]["qc"])
    else:
        started = time.monotonic()
        chromosome_shards, qc = prepare(**prepare_arguments, out=args.out)
        if args.boundary == "cell-batch":
            chromosome_map = json.loads((args.out/"shards.json").read_text())["chromosomes"]
            shared_certificate = (json.loads((args.prepared_reads/"certificate.json").read_text())
                                  if args.prepared_reads is not None else None)
            shards, batch_qc = repartition_by_cb(
                args.out, chromosome_shards, chromosome_map, args.cells_per_batch,
                shared_certificate)
            qc.update(batch_qc)
            (args.out/"shards.json").write_text(json.dumps({
                "boundary": "whole cell within original chromosome",
                "shards": shards, "chromosomes": chromosome_map,
                "batch_to_chromosome_shard": json.loads(
                    (args.out/"batch_chromosomes.json").read_text()
                )["batch_to_chromosome_shard"]}, indent=2, sort_keys=True))
        else:
            shards = chromosome_shards
        outputs = [args.out/"aliases.json", args.out/"bulk.tsv", args.out/"shards.json"]
        if args.boundary == "cell-batch":
            outputs.append(args.out/"batch_chromosomes.json")
            if args.prepared_reads is None:
                outputs += [args.out/"readid_cb.tsv", args.out/"readid_cb.sorted.tsv"]
        outputs += [args.out/shard/name for shard in shards
                    for name in ("annotation.tsv", "reads.bed")]
        write_stage_manifest(prep_manifest, prep_identity, outputs,
                             {"seconds": time.monotonic()-started, "qc": dict(qc),
                              "peak_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss})
    if args.stop_after_prepare:
        marker = args.out/"PREPARED_ONLY.json"
        temporary = marker.with_suffix(".json.tmp")
        temporary.write_text(json.dumps({"complete": True, "scientific_run_complete": False,
            "purpose": "mapping-memory probe before full execution",
            "shards": len(shards), "qc": dict(qc)}, indent=2))
        temporary.replace(marker)
        finish_invocation_ledger(args.out, invocation_id, time.monotonic()-run_started)
        return
    (args.out/"PREPARED_ONLY.json").unlink(missing_ok=True)
    src, commands = args.upstream/"src", []
    if args.boundary == "cell-batch":
        mapping_script, mapping_noop_script = write_mapping_bridges(
            src/"mapping_filtering.R", args.out)
    else:
        mapping_script = src/"mapping_filtering.R"
        mapping_noop_script = None
    shard_payload = json.loads((args.out/"shards.json").read_text())
    batch_to_chromosome = shard_payload.get("batch_to_chromosome_shard", {})
    shared_exons = {}
    parallel_wall_seconds = {}

    def annotate_chromosome(chromosome_shard):
        """Build one complete chromosome exon table for the batch barrier."""
        directory = args.out/chromosome_shard
        identity = stage_identity("chromosome_annotation",
            {"annotation": directory/"annotation.tsv", "bulk": args.out/"bulk.tsv"},
            bound({key: PARAMETERS[key] for key in ("distance_transcriptome", "distance_exon")}))
        manifest = directory/"annotation.complete.json"
        done = completed_stage(manifest, identity) if args.resume else None
        metric = None
        if not done:
            metric = run_command([args.rscript, src/"gtf_processing.R",
                directory/"annotation.tsv", args.out/"bulk.tsv", 600, 30,
                directory/"exons.tsv", directory/"unique.tsv"],
                directory, directory/"annotation.log")
            write_stage_manifest(manifest, identity,
                [directory/"exons.tsv", directory/"unique.tsv"], metric)
        return chromosome_shard, metric

    if args.boundary == "cell-batch":
        chromosome_shards = sorted(set(batch_to_chromosome.values()))
        annotation_results, parallel_wall_seconds["chromosome_annotation"] = bounded_map(
            annotate_chromosome, chromosome_shards, args.r_workers)
        for chromosome_shard, metric in annotation_results:
            if metric is not None:
                commands.append(metric)
            shared_exons[chromosome_shard] = args.out/chromosome_shard/"exons.tsv"

    def map_shard(shard):
        """Run original mapping and IP filtering for one independent complete-CB batch."""
        directory = args.out/shard
        exons_path = (shared_exons[batch_to_chromosome[shard]]
                      if args.boundary == "cell-batch" else directory/"exons.tsv")
        stage_sources = {"annotation": directory/"annotation.tsv",
                         "reads": directory/"reads.bed", "bulk": args.out/"bulk.tsv"}
        if args.boundary == "cell-batch":
            stage_sources["exons"] = exons_path
            stage_sources["mapping_bridge"] = mapping_script
        mapping_stage = "mapping_initial" if args.boundary == "cell-batch" else "preprocess"
        identity = stage_identity(mapping_stage, stage_sources, bound(PARAMETERS))
        manifest = directory/("mapping_initial.complete.json" if args.boundary == "cell-batch"
                             else "preprocess.complete.json")
        done = completed_stage(manifest, identity) if args.resume else None
        filtered_output = directory/("filtered.initial.rds" if args.boundary == "cell-batch"
                                    else "filtered.rds")
        readids_output = directory/("readids.initial.tsv" if args.boundary == "cell-batch"
                                   else "readids.tsv")
        outputs = [filtered_output, readids_output, directory/"diagnostic_provenance.json"]
        if args.boundary == "cell-batch":
            outputs.append(directory/"mapped.rds.concordance_candidates.tsv")
        if args.boundary != "cell-batch":
            outputs.extend([directory/"exons.tsv", directory/"unique.reads"])
        metrics = []
        if not done:
            if args.boundary != "cell-batch":
                metrics.append(run_command([args.rscript, src/"gtf_processing.R",
                    directory/"annotation.tsv", args.out/"bulk.tsv", 600, 30,
                    directory/"exons.tsv", directory/"unique.tsv"], directory,
                    directory/"mapping.log"))
            if (directory/"reads.bed").stat().st_size == 0:
                (directory/"unique.reads").write_text("")
                readids_output.write_text("")
                filtered_output.write_text("NO_EVIDENCE\n")
                install_diagnostic_sink(directory)
                metrics.append({"status": "no_evidence", "seconds": 0,
                                "peak_rss_kib": 0})
            else:
                (directory/"empty_ip.tsv").write_text(
                    "seqnames.ip\tstart.ip\tend.ip\tc4\tc5\tstrand.ip\n")
                install_diagnostic_sink(directory)
                metrics.append(run_command([args.rscript, mapping_script,
                    directory/"reads.bed", exons_path, directory/"mapped.rds"],
                    directory, directory/"mapping.log"))
                local_unique = (Path("/dev/null") if args.boundary == "cell-batch"
                                else directory/"unique.reads")
                metrics.append(run_command([args.rscript, src/"ip_filtering.R",
                    directory/"mapped.rds", directory/"empty_ip.tsv", 60,
                    readids_output, local_unique,
                    directory/"mapped.ip", filtered_output],
                    directory, directory/"mapping.log"))
            if not filtered_output.is_file() or filtered_output.stat().st_size == 0:
                raise RuntimeError(f"Original filtering did not publish a usable RDS: {directory}")
            if not readids_output.is_file():
                raise RuntimeError(f"Original filtering did not publish read IDs: {directory}")
            write_stage_manifest(manifest, identity, outputs, {"commands": metrics})
            (directory/"mapped.rds").unlink(missing_ok=True)
        return shard, metrics

    mapping_results, parallel_wall_seconds["mapping_initial"] = bounded_map(
        map_shard, shards, args.r_workers)
    for shard, metrics in mapping_results:
        commands.extend(metrics)
    unique_files = ([args.out/shard/"unique.reads" for shard in shards]
                    if args.boundary != "cell-batch" else [])
    if args.boundary == "cell-batch":
        by_chromosome = defaultdict(list)
        for shard in shards:
            by_chromosome[batch_to_chromosome[shard]].append(shard)
        for chromosome_shard, chromosome_batches in sorted(by_chromosome.items()):
            chromosome_directory = args.out/chromosome_shard
            candidate_paths = [args.out/shard/"mapped.rds.concordance_candidates.tsv"
                               for shard in chromosome_batches]
            correction_sources = {f"candidates_{index:05d}": path
                                  for index, path in enumerate(candidate_paths)}
            correction_sources["noop_bridge"] = mapping_noop_script
            identity = stage_identity("chromosome_concordance_bug", correction_sources,
                bound({"semantics": "native one-column data.table membership: singleton only"}))
            manifest = chromosome_directory/"concordance.complete.json"
            done = completed_stage(manifest, identity) if args.resume else None
            global_candidates = chromosome_directory/"concordance_candidates.tsv"
            provenance = chromosome_directory/"concordance_provenance.json"
            if not done:
                plan = concordance_correction_plan(candidate_paths, global_candidates, args.out)
                metrics = []
                rerun_indexes = set(plan["rerun_indexes"])
                for index, shard in enumerate(chromosome_batches):
                    if index in rerun_indexes:
                        continue
                    directory = args.out/shard
                    for initial_name, final_name in (("filtered.initial.rds", "filtered.rds"),
                                                     ("readids.initial.tsv", "readids.tsv")):
                        temporary = directory/(final_name + ".tmp")
                        temporary.unlink(missing_ok=True)
                        temporary.hardlink_to(directory/initial_name)
                        temporary.replace(directory/final_name)
                for index in plan["rerun_indexes"]:
                    directory = args.out/chromosome_batches[index]
                    mapped = directory/"mapped.corrected.rds"
                    readids = directory/"readids.corrected.tsv"
                    filtered = directory/"filtered.corrected.rds"
                    metrics.append(run_command([args.rscript, mapping_noop_script,
                        directory/"reads.bed", shared_exons[chromosome_shard], mapped],
                        directory, args.out/"run.log"))
                    metrics.append(run_command([args.rscript, src/"ip_filtering.R", mapped,
                        directory/"empty_ip.tsv", 60, readids, Path("/dev/null"),
                        directory/"mapped.corrected.ip", filtered],
                        directory, args.out/"run.log"))
                    if not filtered.is_file() or filtered.stat().st_size == 0 or not readids.is_file():
                        raise RuntimeError(f"Concordance correction failed: {directory}")
                    filtered.replace(directory/"filtered.rds")
                    readids.replace(directory/"readids.tsv")
                    mapped.unlink(missing_ok=True)
                provenance.write_text(json.dumps({
                    "upstream_bug": "ftrs %in% one-column data.table matches only one-row tables",
                    "preserved_semantics": "delete singleton iff chromosome-global distinct count is one",
                    **plan}, indent=2, sort_keys=True))
                outputs = [global_candidates, provenance]
                outputs += [args.out/shard/name for shard in chromosome_batches
                            for name in ("filtered.rds", "readids.tsv")]
                write_stage_manifest(manifest, identity, outputs, {"commands": metrics, "qc": plan})
                commands.extend(metrics)
        pair_files = []
        for shard in shards:
            directory = args.out/shard
            identity = stage_identity("observed_pairs", {"filtered": directory/"filtered.rds"},
                                      bound({"scope": "post-original-filtering"}))
            manifest = directory/"observed_pairs.complete.json"
            done = completed_stage(manifest, identity) if args.resume else None
            if not done:
                metric = run_command([args.rscript, Path(__file__).with_name("batch_pairs.R"),
                    directory/"filtered.rds", directory/"observed_pairs.tsv"],
                    directory, args.out/"run.log")
                write_stage_manifest(manifest, identity, [directory/"observed_pairs.tsv"], metric)
                commands.append(metric)
            pair_files.append(directory/"observed_pairs.tsv")
        identity = stage_identity("global_observed_pairs",
            {f"batch_{index:05d}": path for index, path in enumerate(pair_files)},
            bound({"definition": "observed global uniqueN(transcript_name)==1"}))
        manifest = args.out/"global_observed_pairs.complete.json"
        done = completed_stage(manifest, identity) if args.resume else None
        if not done:
            started = time.monotonic()
            pair_qc = reconcile_observed_pairs(
                pair_files, args.out/"global_unique_genes.tsv", args.out)
            metric = {"seconds": time.monotonic()-started, "peak_rss_kib":
                      resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, "qc": pair_qc}
            write_stage_manifest(manifest, identity,
                [args.out/"observed_pairs.tsv", args.out/"global_unique_genes.tsv"], metric)
        for shard in shards:
            directory = args.out/shard
            identity = stage_identity("global_unique_rows",
                {"filtered": directory/"filtered.rds",
                 "global_unique_genes": args.out/"global_unique_genes.tsv"},
                bound({"columns": "original ip_filtering unique.reads"}))
            manifest = directory/"global_unique_rows.complete.json"
            done = completed_stage(manifest, identity) if args.resume else None
            if not done:
                metric = run_command([args.rscript, Path(__file__).with_name("batch_unique.R"),
                    directory/"filtered.rds", args.out/"global_unique_genes.tsv",
                    directory/"unique.reads"], directory, args.out/"run.log")
                write_stage_manifest(manifest, identity, [directory/"unique.reads"], metric)
                commands.append(metric)
            unique_files.append(directory/"unique.reads")
    concatenate(unique_files, args.out/"all_unique.reads")
    if args.calibration_probabilities is None:
        cal_sources = {"unique_reads": args.out/"all_unique.reads"}
    else:
        cal_sources = {"provided_probabilities": args.calibration_probabilities}
    cal_identity = stage_identity("calibration", cal_sources,
        bound({key: PARAMETERS[key] for key in ("gene_fraction","binsize")}))
    cal_manifest = args.out/"calibration.complete.json"
    done = completed_stage(cal_manifest, cal_identity) if args.resume else None
    if not done:
        if args.calibration_probabilities is None:
            metric = run_command([args.rscript, src/"compute_prob.R",
                args.out/"all_unique.reads", "98%", 20, args.out/"probabilities.tsv",
                args.out/"probabilities.pdf"], args.out, args.out/"run.log")
            outputs = [args.out/"probabilities.tsv", args.out/"probabilities.pdf"]
        else:
            shutil.copyfile(args.calibration_probabilities, args.out/"probabilities.tsv")
            metric = {"status": "reference_calibrated", "seconds": 0,
                      "peak_rss_kib": 0}
            outputs = [args.out/"probabilities.tsv"]
            (args.out/"calibration_provenance.json").write_text(json.dumps(
                {"source": file_fingerprint(args.calibration_probabilities),
                 "mode": "reference-calibrated capture distribution"}, indent=2))
            outputs.append(args.out/"calibration_provenance.json")
        write_stage_manifest(cal_manifest, cal_identity, outputs, metric)
        commands.append(metric)
    def fragment_shard(shard):
        """Compute original fragment weights for one finalized independent batch."""
        directory = args.out/shard
        no_filtered_evidence = (args.boundary == "cell-batch" and not
                                filtered_batch_has_evidence(directory/"observed_pairs.tsv"))
        if (directory/"reads.bed").stat().st_size == 0 or no_filtered_evidence:
            (directory/"no_filtered_evidence.json").write_text(json.dumps({
                "status": "zero rows after successful original filtering",
                "annotation_retained": True}, indent=2))
            return shard, None
        identity = stage_identity("fragment_probabilities",
            {"filtered": directory/"filtered.rds",
             "probabilities": args.out/"probabilities.tsv"}, bound(PARAMETERS))
        manifest = directory/"fragments.complete.json"
        done = completed_stage(manifest, identity) if args.resume else None
        metric = None
        if not done:
            metric = run_command([args.rscript, src/"fragment_probabilities.R",
                directory/"filtered.rds", args.out/"probabilities.tsv",
                directory/"fragments.tsv"], directory, directory/"fragments.log")
            write_stage_manifest(manifest, identity,
                                 [directory/"fragments.tsv"], metric)
        return shard, metric

    fragment_results, parallel_wall_seconds["fragment_probabilities"] = bounded_map(
        fragment_shard, shards, args.r_workers)
    fragment_shards = []
    for shard, metric in fragment_results:
        if metric is not None:
            commands.append(metric)
        if (args.out/shard/"fragments.tsv").is_file():
            fragment_shards.append(shard)
    fragment_files = [args.out/shard/"fragments.tsv" for shard in fragment_shards]

    def em_shard(shard):
        """Run the unchanged original parsed EM for one weighted independent batch."""
        directory = args.out/shard
        identity = stage_identity("em",
            {"fragments": directory/"fragments.tsv",
             "em_source": args.upstream/"src/em_algorithm.R"},
            bound({"MAX_IT": 30, "round_digits": 3}))
        manifest = directory/"em.complete.json"
        done = completed_stage(manifest, identity) if args.resume else None
        metric = None
        if not done:
            metric = run_command([args.rscript, Path(__file__).with_name("em_batch.R"),
                args.upstream, directory/"fragments.tsv",
                directory/"estimates_aliased.tsv", directory/"support_aliased.tsv"],
                directory, directory/"em.log")
            write_stage_manifest(manifest, identity,
                [directory/"estimates_aliased.tsv", directory/"support_aliased.tsv"], metric)
        return shard, metric

    em_results, parallel_wall_seconds["em"] = bounded_map(
        em_shard, fragment_shards, args.r_workers)
    for shard, metric in em_results:
        if metric is not None:
            commands.append(metric)
    estimate_files = [args.out/shard/"estimates_aliased.tsv" for shard in fragment_shards]
    support_files = [args.out/shard/"support_aliased.tsv" for shard in fragment_shards]
    aliases = json.loads((args.out/"aliases.json").read_text())
    for name, paths in (("estimates", estimate_files), ("support", support_files)):
        destination = args.out/f"{name}.tsv"
        temporary = destination.with_suffix(destination.suffix+".tmp")
        wrote_header = False
        with temporary.open("w", newline="") as target:
            writer = None
            for path in paths:
                with path.open() as source:
                    reader = csv.DictReader(source, delimiter="\t")
                    if reader.fieldnames is None:
                        continue
                    if writer is None:
                        writer = csv.DictWriter(target, fieldnames=reader.fieldnames,
                                                delimiter="\t", lineterminator="\n")
                        writer.writeheader(); wrote_header = True
                    for row in reader:
                        for key in ("gene_name", "transcript_name"):
                            if key in row:
                                row[key] = aliases[row[key]]
                        writer.writerow(row)
        if not wrote_header:
            raise RuntimeError(f"No completed nonempty SCALPEL {name} shards")
        temporary.replace(destination)
    records = collect_stage_records(args.out)
    resources = summarize_resources(records)
    (args.out/"stage_metrics.json").write_text(json.dumps(records, indent=2))
    (args.out/"commands.json").write_text(json.dumps(resources["commands"], indent=2))
    adoption_path = args.out/"prepare_adoption.json"
    adopted_wall_seconds = (float(json.loads(adoption_path.read_text())
                                  ["adopted_wall_seconds"])
                            if adoption_path.is_file() else 0.0)
    invocation_wall_seconds = time.monotonic()-run_started
    finish_invocation_ledger(args.out, invocation_id, invocation_wall_seconds)
    summary = {"upstream_commit": revision, "license": "AGPL-3.0",
        "inputs": {key: str(value) for key, value in vars(args).items()},
        "parameters": PARAMETERS, "qc": dict(qc), "shards": len(shards),
        "completed_fragment_shards": len(fragment_files),
        "wall_seconds": (adopted_wall_seconds + prior_invocation_seconds
                         + invocation_wall_seconds),
        "wall_seconds_complete": incomplete_invocations == 0,
        "unclosed_prior_invocations": incomplete_invocations,
        "runtime_scope": ("validated adopted preparation plus completed invocations; "
                          "unclosed interrupted attempts have unknown duration and are "
                          "flagged instead of silently estimated"),
        "adopted_wall_seconds": adopted_wall_seconds,
        "summed_command_wall_seconds": resources["stage_wall_seconds"],
        "parallel_stage_wall_seconds": parallel_wall_seconds,
        "wall_seconds_this_invocation": invocation_wall_seconds,
        "peak_child_rss_kib": resources["peak_child_rss_kib"],
        "internal_priming": "disabled, original empty-IP branch",
        "execution_boundary": args.boundary,
        "r_concurrency": args.r_workers,
        "convergence": "Upstream MAX_IT=30, tolerance=0.01; iteration counts not exposed",
        "count_scaling": "rounded upstream relative abundance times positive-probability cell/gene molecules",
        "resume": "stages bound to SHA-256 inputs, parameters, harness, and upstream sources"}
    temporary = args.out/"summary.json.tmp"
    temporary.write_text(json.dumps(summary, indent=2))
    temporary.replace(args.out/"summary.json")


if __name__ == "__main__":
    main()
