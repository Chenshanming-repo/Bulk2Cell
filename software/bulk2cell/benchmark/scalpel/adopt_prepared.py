#!/usr/bin/env python3
"""Explicitly adopt a verified prepared SCALPEL stage into a reviewed runner."""
import argparse
import ast
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import time

PREPARATION_SYMBOLS = (
    "HandlePool", "assign_cb_batches", "_cell_from_fragment", "upstream_read_id",
    "repartition_by_cb", "_blocks", "prepare", "validate_prepared_reads")


def sha256(path):
    """Return the streaming SHA-256 digest of a file."""
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def lineage_matches_fingerprint(lineage, fingerprint):
    """Compare lineage bytes with a runner fingerprint's equivalent size field."""
    return (lineage.get("path") == fingerprint.get("path") and
            lineage.get("sha256") == fingerprint.get("sha256") and
            lineage.get("bytes") == fingerprint.get("size"))


def verified_lineage_chain(lineage_path, fingerprint, max_depth=16):
    """Verify every certificate hop and return fingerprints of the exact chain."""
    path = Path(lineage_path).resolve()
    visited, evidence = set(), []
    for _ in range(max_depth):
        if path in visited:
            raise ValueError("Prepared-read lineage contains a cycle")
        visited.add(path)
        evidence.append(sha256_record(path))
        lineage = json.loads(path.read_text())
        source_bundle = lineage.get("source_bundle")
        if not source_bundle:
            return None
        certificate_path = Path(source_bundle).resolve()/"certificate.json"
        if not certificate_path.is_file():
            return None
        actual = sha256_record(certificate_path)
        declared = lineage.get("source_certificate", {})
        if not lineage_matches_fingerprint(declared, actual):
            raise ValueError("Prepared-read lineage certificate hop is invalid")
        evidence.append(actual)
        if (actual["path"] == fingerprint.get("path") and
                actual["sha256"] == fingerprint.get("sha256") and
                actual["size"] == fingerprint.get("size")):
            return evidence
        path = Path(source_bundle).resolve()/"lineage.json"
        if not path.is_file():
            return None
    raise ValueError("Prepared-read lineage exceeds the maximum depth")


def lineage_contains_fingerprint(lineage_path, fingerprint, max_depth=16):
    """Return whether a fully verified lineage reaches an exact fingerprint."""
    return verified_lineage_chain(lineage_path, fingerprint, max_depth) is not None


def sha256_record(path):
    """Return path, size, and streaming SHA-256 for lineage evidence."""
    path = Path(path).resolve()
    return {"path": str(path), "size": path.stat().st_size, "sha256": sha256(path)}


def load_runner(path, name):
    """Load a runner module from an explicitly selected source path."""
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def symbol_hashes(path):
    """Hash exact source segments for preparation functions and classes."""
    source = Path(path).read_text()
    tree = ast.parse(source)
    nodes = {node.name: node for node in tree.body
             if isinstance(node, (ast.FunctionDef, ast.ClassDef))}
    missing = set(PREPARATION_SYMBOLS) - set(nodes)
    if missing:
        raise ValueError(f"Missing preparation symbols: {sorted(missing)}")
    return {name: hashlib.sha256(
        ast.get_source_segment(source, nodes[name]).encode()).hexdigest()
        for name in PREPARATION_SYMBOLS}


def run_binding(runner, upstream, revision):
    """Build the exact run identity used by the candidate runner."""
    directory = Path(runner.__file__).parent
    return {
        "runner": runner.file_fingerprint(runner.__file__),
        "em_batch": runner.file_fingerprint(directory/"em_batch.R"),
        "batch_pairs": runner.file_fingerprint(directory/"batch_pairs.R"),
        "batch_unique": runner.file_fingerprint(directory/"batch_unique.R"),
        "upstream_manifest": runner.file_fingerprint(directory/"upstream-sha256.json"),
        "upstream_commit": revision}


def main():
    """Validate old outputs and atomically publish an explicit adoption manifest."""
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ("out", "old-run", "candidate-run", "upstream", "gtf", "quant",
                 "bam", "barcodes", "prepared-reads"):
        parser.add_argument("--"+name, required=True, type=Path)
    parser.add_argument("--cells-per-batch", required=True, type=int)
    args = parser.parse_args()
    for key, value in vars(args).items():
        if isinstance(value, Path):
            setattr(args, key, value.resolve())
    if args.cells_per_batch < 1:
        parser.error("--cells-per-batch must be at least 1")
    adoption_path = args.out/"prepare_adoption.json"
    archive_path = args.out/"prepare.pre_workers.archived.json"
    completion_path = args.out/"PREPARE_ADOPTION_COMPLETE.json"
    manifest_path = args.out/"prepare.complete.json"
    prepared_only = json.loads((args.out/"PREPARED_ONLY.json").read_text())
    if not prepared_only.get("complete") or prepared_only.get("scientific_run_complete"):
        raise ValueError("Old preparation completion marker is invalid")
    candidate = load_runner(args.candidate_run, "scalpel_candidate")
    old_runner = load_runner(args.old_run, "scalpel_old")
    lock_handle = candidate.acquire_run_lock(args.out)
    legacy_path = archive_path if archive_path.is_file() else manifest_path
    old_manifest = json.loads(legacy_path.read_text())
    old_identity = {key: old_manifest[key]
                    for key in ("stage", "sources", "parameters")}
    if not old_runner.completed_stage(legacy_path, old_identity):
        raise ValueError("Old prepare manifest or one of its declared outputs is invalid")
    declared_runner = old_manifest["parameters"].get("run_binding", {}).get("runner")
    supplied_runner = candidate.file_fingerprint(args.old_run)
    if declared_runner is None or any(declared_runner.get(key) != supplied_runner[key]
                                      for key in ("size", "sha256")):
        raise ValueError("Supplied old runner did not produce the legacy manifest")
    if old_manifest["parameters"].get("boundary") != "cell-batch":
        raise ValueError("Old preparation boundary is not complete-CB batching")
    if old_manifest["parameters"].get("cells_per_batch") != args.cells_per_batch:
        raise ValueError("Requested cell batch size differs from old preparation")
    requested = {"gtf": args.gtf, "quant": args.quant, "bam": args.bam,
                 "barcodes": args.barcodes}
    for name, path in requested.items():
        if old_manifest["sources"].get(name) != candidate.file_fingerprint(path):
            raise ValueError(f"Requested {name} differs from old prepared input")
    lineage_path = args.prepared_reads/"lineage.json"
    old_certificate = old_manifest["sources"].get("prepared_certificate")
    if old_certificate is None:
        raise ValueError("Old preparation was not built from certified shared reads")
    lineage_evidence = verified_lineage_chain(lineage_path, old_certificate)
    if lineage_evidence is None:
        raise ValueError("Versioned prepared-read lineage differs from old certificate")
    old_hashes, new_hashes = symbol_hashes(args.old_run), symbol_hashes(args.candidate_run)
    if old_hashes != new_hashes:
        changed = sorted(name for name in old_hashes if old_hashes[name] != new_hashes[name])
        raise ValueError(f"Preparation implementation changed: {changed}")
    shards = json.loads((args.out/"shards.json").read_text())
    if shards.get("boundary") != "whole cell within original chromosome":
        raise ValueError("Only complete-CB preparation may be adopted")
    chromosomes = set(shards["chromosomes"].values())
    prepared = candidate.validate_prepared_reads(
        args.prepared_reads, args.bam, args.barcodes, chromosomes)
    if prepared["certificate"]["cross_cb_read_ids"] != 0:
        raise ValueError("Versioned certificate has cross-CB read IDs")
    revision = candidate.verify_upstream(args.upstream)
    binding = run_binding(candidate, args.upstream, revision)
    inputs = {**requested,
              "prepared_certificate": args.prepared_reads/"certificate.json"}
    identity = candidate.stage_identity("prepare", inputs, {
        "boundary": "cell-batch", "cells_per_batch": args.cells_per_batch,
        "run_binding": binding})
    old_outputs = [Path(row) for row in old_manifest["outputs"]]
    if any(args.out not in path.parents for path in old_outputs):
        raise ValueError("Old prepare manifest declares output outside run directory")
    if not archive_path.exists():
        temporary = archive_path.with_suffix(".json.tmp")
        shutil.copy2(manifest_path, temporary)
        temporary.replace(archive_path)
    adoption = {
        "schema": "bulk2cell.scalpel.prepare_adoption.v1",
        "status": "ready_for_manifest",
        "old_manifest": candidate.file_fingerprint(archive_path),
        "old_runner": candidate.file_fingerprint(args.old_run),
        "candidate_runner": candidate.file_fingerprint(args.candidate_run),
        "adopter": candidate.file_fingerprint(__file__),
        "preparation_symbol_sha256": new_hashes,
        "new_prepared_certificate": candidate.file_fingerprint(
            args.prepared_reads/"certificate.json"),
        "versioned_lineage": candidate.file_fingerprint(lineage_path),
        "verified_lineage_chain": lineage_evidence,
        "declared_outputs_revalidated": len(old_outputs),
        "adopted_output_fingerprints": [candidate.file_fingerprint(path)
                                        for path in old_outputs],
        "cells_per_batch": args.cells_per_batch,
        "shards": len(shards["shards"]),
        "adopted_wall_seconds": float(old_manifest["metrics"]["seconds"]),
        "adopted_unix": time.time(),
        "reason": "reviewed bounded R workers; preparation implementation byte-identical"
    }
    temporary = adoption_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(adoption, indent=2)+"\n")
    temporary.replace(adoption_path)
    candidate.write_stage_manifest(
        manifest_path, identity, old_outputs+[archive_path, adoption_path],
        {**old_manifest["metrics"], "adoption": adoption})
    completion = {"completed": True,
                  "adoption": candidate.file_fingerprint(adoption_path),
                  "prepare_manifest": candidate.file_fingerprint(manifest_path)}
    temporary = completion_path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(completion, indent=2)+"\n")
    temporary.replace(completion_path)
    if not candidate.completed_stage(manifest_path, identity):
        raise RuntimeError("New prepare manifest failed post-publication verification")
    candidate.validate_prepare_adoption(args.out)
    print(json.dumps({"adopted": True, "outputs": len(old_outputs),
                      "manifest": str(manifest_path)}, indent=2))


if __name__ == "__main__":
    main()
