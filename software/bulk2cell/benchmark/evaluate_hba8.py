#!/usr/bin/env python3
"""Evaluate completed, barcode-matched bulk2cell and SCALPEL runs.

The script treats Cell Ranger clusters as a reference partition, never as
isoform-level truth.  It validates method-specific completion records, loads
cells-by-isoform sparse matrices with explicit identifiers, reports native
callability and resource use, and evaluates every method on the intersection
of active barcodes.  Outputs are strict JSON, TSV, and standalone PNG figures.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import textwrap
from typing import Mapping, NamedTuple, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse

from bulk2cell.evaluation import (
    evaluate_cells,
    load_matrix_market,
    read_cell_ranger_clusters,
)



class BenchmarkRun(NamedTuple):
    """One validated completed run and its explicitly ordered sparse matrix."""

    label: str
    path: Path
    kind: str
    matrix: sparse.csr_matrix
    barcodes: tuple[str, ...]
    transcript_ids: tuple[str, ...]
    gene_ids: tuple[str, ...]
    manifest: Mapping[str, object]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    """Read a headered TSV, rejecting a missing header."""

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"TSV has no header: {path}")
        return list(reader)


def _write_tsv(path: Path, header: Sequence[str], rows: Sequence[Sequence[object]]) -> None:
    """Write a deterministic headered TSV."""

    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def load_run(label: str, path: str | Path, kind: str) -> BenchmarkRun:
    """Load one completed run after validating its method-specific manifest."""

    directory = Path(path)
    if kind not in {"bulk2cell", "scalpel"}:
        raise ValueError("kind must be 'bulk2cell' or 'scalpel'")
    if kind == "bulk2cell":
        manifest_path = directory / "run.json"
        manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
        if manifest.get("status") != "complete":
            raise ValueError(f"{directory} is not a completed bulk2cell run")
    else:
        manifest_path = directory / "summary.json"
        matrix_summary = directory / "matrix_summary.json"
        if not manifest_path.exists() or not matrix_summary.exists():
            raise ValueError(f"{directory} is not a completed exported SCALPEL run")
        manifest = json.loads(manifest_path.read_text())
    barcode_rows = _read_tsv(directory / "barcodes.tsv")
    isoform_rows = _read_tsv(directory / "isoforms.tsv")
    barcodes = tuple(row["barcode"] for row in barcode_rows)
    transcripts = tuple(row["transcript_id"] for row in isoform_rows)
    genes = tuple(row["gene_id"] for row in isoform_rows)
    matrix, _, _ = load_matrix_market(
        directory / "isoform_counts.mtx.gz", barcodes, transcripts
    )
    return BenchmarkRun(label, directory, kind, matrix, barcodes, transcripts, genes, manifest)


def summarize_run(run: BenchmarkRun) -> tuple[dict[str, object], list[tuple[object, ...]]]:
    """Compute native-cell callability, resources, convergence, and gene totals."""

    cell_totals = np.asarray(run.matrix.sum(axis=1)).ravel()
    unique_genes = tuple(sorted(set(run.gene_ids)))
    gene_index = {gene: index for index, gene in enumerate(unique_genes)}
    columns = np.asarray([gene_index[gene] for gene in run.gene_ids])
    projector = sparse.csr_matrix(
        (np.ones(len(columns)), (np.arange(len(columns)), columns)),
        shape=(len(columns), len(unique_genes)),
    )
    gene_matrix = run.matrix @ projector
    gene_totals = np.asarray(gene_matrix.sum(axis=0)).ravel()
    gene_rows = [(run.label, gene, float(total)) for gene, total in zip(unique_genes, gene_totals)]
    if run.kind == "bulk2cell":
        qc = run.manifest.get("qc", {})
        if not isinstance(qc, Mapping):
            qc = {}
        convergence = {
            "status": "reported_by_bulk2cell",
            "nonconverged_cell_gene_fits": int(qc.get("em_nonconverged_cells", 0)),
        }
        runtime = float(run.manifest["elapsed_seconds"])
        peak_rss = float(run.manifest["peak_rss_mib"])
    else:
        convergence = {
            "status": "not_exposed_by_upstream_scalpel",
            "detail": run.manifest.get("convergence"),
        }
        runtime = float(run.manifest["wall_seconds"])
        peak_rss = float(run.manifest["peak_child_rss_kib"]) / 1024.0
    summary: dict[str, object] = {
        "method": run.label,
        "kind": run.kind,
        "whitelist_cells": len(run.barcodes),
        "active_cells": int(np.count_nonzero(cell_totals > 0)),
        "inactive_cells": int(np.count_nonzero(cell_totals == 0)),
        "isoform_features": len(run.transcript_ids),
        "genes": len(unique_genes),
        "callable_cell_gene_pairs": int(gene_matrix.getnnz(axis=1).sum()),
        "total_estimated_counts": float(run.matrix.sum()),
        "runtime_seconds": runtime,
        "peak_rss_mib": peak_rss,
        "convergence": convergence,
        "runtime_scope": run.manifest.get("runtime_scope", "independent_method_run"),
        "peak_rss_scope": run.manifest.get("peak_rss_scope",
            "parent_only" if run.kind == "bulk2cell" else "child_lifetime_maximum_not_aggregate"),
        "external_training_elapsed_seconds": run.manifest.get("external_training_elapsed_seconds", 0.0),
        "training_metrics": run.manifest.get("training_metrics", {}),
        "fit_metrics": run.manifest.get("fit_metrics", {}),
    }
    return summary, gene_rows


def validate_matched_whitelists(runs: Sequence[BenchmarkRun]) -> None:
    """Require exact whitelist set equality before shared-active filtering."""

    if not runs:
        raise ValueError("at least one run is required")
    expected = set(runs[0].barcodes)
    for run in runs[1:]:
        if set(run.barcodes) != expected:
            raise ValueError(f"whitelist barcode sets differ for {run.label}")


def shared_active_barcodes(runs: Sequence[BenchmarkRun]) -> tuple[str, ...]:
    """Return sorted barcodes with positive counts in every supplied run."""

    if not runs:
        raise ValueError("at least one run is required")
    active_sets = []
    for run in runs:
        totals = np.asarray(run.matrix.sum(axis=1)).ravel()
        active_sets.append({barcode for barcode, total in zip(run.barcodes, totals) if total > 0})
    return tuple(sorted(set.intersection(*active_sets)))


def _safe_label(label: str) -> str:
    """Convert a method label to a stable filename stem."""

    safe = "".join(character if character.isalnum() else "_" for character in label)
    return safe.strip("_").lower()


def _scatter(path: Path, coordinates: np.ndarray, labels: Sequence[object], title: str,
             axes: tuple[str, str]) -> None:
    """Write a standalone reference-colored two-dimensional embedding plot."""

    categories = sorted(set(map(str, labels)), key=lambda value: (0, int(value)) if value.isdigit() else (1, value))
    encoded = np.asarray([categories.index(str(value)) for value in labels])
    figure, axis = plt.subplots(figsize=(7.2, 6.0), constrained_layout=True)
    points = axis.scatter(coordinates[:, 0], coordinates[:, 1], c=encoded, s=5,
                          cmap="tab20", alpha=0.75, linewidths=0, rasterized=True)
    axis.set_xlabel(axes[0])
    axis.set_ylabel(axes[1])
    axis.set_title(textwrap.fill(title, width=60))
    handles, _ = points.legend_elements(num=None)
    if len(categories) <= 20:
        axis.legend(handles, categories, title="Cell Ranger cluster", fontsize=7,
                    loc="best", frameon=False)
    figure.savefig(path, dpi=220)
    plt.close(figure)


def evaluate_shared(
    run: BenchmarkRun,
    shared_barcodes: Sequence[str],
    reference_labels: Mapping[str, object],
    output: str | Path,
    *,
    n_clusters: int,
    seed: int,
    compute_umap: bool,
) -> dict[str, object]:
    """Evaluate one run on a fixed shared barcode set and write coordinates/plots."""

    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    positions = {barcode: index for index, barcode in enumerate(run.barcodes)}
    missing = [barcode for barcode in shared_barcodes if barcode not in positions]
    if missing:
        raise ValueError(f"{run.label} lacks {len(missing)} shared barcodes")
    matrix = run.matrix[[positions[barcode] for barcode in shared_barcodes]]
    result = evaluate_cells(
        matrix,
        shared_barcodes,
        n_clusters,
        reference_labels=reference_labels,
        random_state=seed,
        compute_umap=compute_umap,
    )
    stem = _safe_label(run.label)
    _write_tsv(
        output / f"{stem}_shared_clusters.tsv",
        ["barcode", "inferred_cluster", "cell_ranger_cluster"],
        [(barcode, int(cluster), reference_labels[barcode])
         for barcode, cluster in zip(shared_barcodes, result.clusters)],
    )
    _write_tsv(
        output / f"{stem}_shared_embedding.tsv",
        ["barcode", *[f"SVD{i + 1}" for i in range(result.embedding.shape[1])]],
        [(barcode, *row) for barcode, row in zip(shared_barcodes, result.embedding)],
    )
    colors = [reference_labels[barcode] for barcode in shared_barcodes]
    _scatter(output / f"{stem}_shared_svd.png", result.embedding[:, :2], colors,
             f"{run.label}: shared active cells in SVD space", ("SVD1", "SVD2"))
    if result.umap is not None:
        _write_tsv(
            output / f"{stem}_shared_umap.tsv", ["barcode", "UMAP1", "UMAP2"],
            [(barcode, *row) for barcode, row in zip(shared_barcodes, result.umap)],
        )
        _scatter(output / f"{stem}_shared_umap.png", result.umap, colors,
                 f"{run.label}: shared active cells", ("UMAP1", "UMAP2"))
    return {
        "method": run.label,
        "shared_active_cells": len(shared_barcodes),
        "n_clusters": n_clusters,
        "svd_components": int(result.embedding.shape[1]),
        "ari_cell_ranger_reference": result.adjusted_rand,
        "nmi_cell_ranger_reference": result.normalized_mutual_info,
        "silhouette_svd": result.silhouette,
    }


def _parse_run(value: str) -> tuple[str, str, Path]:
    """Parse ``LABEL:KIND:PATH`` while allowing colons inside the path."""

    fields = value.split(":", 2)
    if len(fields) != 3:
        raise argparse.ArgumentTypeError("run must be LABEL:KIND:PATH")
    return fields[0], fields[1], Path(fields[2])


def main(argv: Sequence[str] | None = None) -> int:
    """Run the matched benchmark evaluation from command-line arguments."""

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", action="append", required=True, type=_parse_run,
                        metavar="LABEL:KIND:PATH")
    parser.add_argument("--clusters", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--n-clusters", type=int, default=10)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--no-umap", action="store_true")
    parser.add_argument("--scope", choices=("panel", "full-sample"), default="panel",
                        help="Declared input scope; gene counts are derived from feature metadata")
    args = parser.parse_args(argv)
    if args.out.exists():
        raise FileExistsError(f"output already exists: {args.out}")
    runs = [load_run(label, path, kind) for label, kind, path in args.run]
    validate_matched_whitelists(runs)
    labels = read_cell_ranger_clusters(args.clusters)
    shared = shared_active_barcodes(runs)
    if len(shared) <= args.n_clusters:
        raise ValueError("shared active cell count must exceed n-clusters")
    args.out.mkdir(parents=True)
    summaries, gene_rows, metrics = [], [], []
    for run in runs:
        summary, rows = summarize_run(run)
        summaries.append(summary)
        gene_rows.extend(rows)
        metrics.append(evaluate_shared(
            run, shared, labels, args.out, n_clusters=args.n_clusters,
            seed=args.seed, compute_umap=not args.no_umap,
        ))
    _write_tsv(
        args.out / "run_summary.tsv",
        ["method", "kind", "whitelist_cells", "active_cells", "inactive_cells",
         "isoform_features", "genes", "callable_cell_gene_pairs",
         "total_estimated_counts", "runtime_seconds", "peak_rss_mib", "runtime_scope", "peak_rss_scope"],
        [[row[key] for key in (
            "method", "kind", "whitelist_cells", "active_cells", "inactive_cells",
            "isoform_features", "genes", "callable_cell_gene_pairs",
            "total_estimated_counts", "runtime_seconds", "peak_rss_mib", "runtime_scope", "peak_rss_scope",
        )] for row in summaries],
    )
    _write_tsv(args.out / "gene_totals.tsv", ["method", "gene_id", "estimated_count"], gene_rows)
    _write_tsv(
        args.out / "shared_clustering_metrics.tsv",
        ["method", "shared_active_cells", "n_clusters", "svd_components",
         "ari_cell_ranger_reference", "nmi_cell_ranger_reference", "silhouette_svd"],
        [[row[key] for key in (
            "method", "shared_active_cells", "n_clusters", "svd_components",
            "ari_cell_ranger_reference", "nmi_cell_ranger_reference", "silhouette_svd",
        )] for row in metrics],
    )
    gene_sets = [set(run.gene_ids) for run in runs]
    genes_union = len(set.union(*gene_sets))
    genes_shared = len(set.intersection(*gene_sets))
    scope_limits = ([f"The targeted {genes_union}-gene panel limits clustering and generalization."]
                    if args.scope == "panel" else
                    ["Full-sample scope retains catalog eligibility and method-specific feature/filter limitations."])
    report = {
        "status": "complete",
        "whitelist_validation": {"status": "exact_set_equality", "cells_per_run": len(runs[0].barcodes)},
        "shared_active_cells": len(shared),
        "shared_barcode_definition": "positive total estimated count in every evaluated method",
        "cell_ranger_role": "reference clustering, not biological or isoform truth",
        "independent_truth": False,
        "scope": args.scope,
        "panel_genes": genes_union,
        "genes_union": genes_union,
        "genes_shared": genes_shared,
        "limitations": [
            "No independent biological truth or held-out read split was used.",
            "Cell Ranger clusters are a reference partition, not truth.",
            "Long-read annotation, priors, and TES evidence are not independent validation.",
            "SCALPEL calibration and output-scaling constraints require method-specific interpretation.",
            "Joint bulk2cell variant wall times are shared; reported process RSS maxima are not aggregate memory.",
        ] + scope_limits,
        "runs": summaries,
        "shared_clustering": metrics,
    }
    (args.out / "evaluation.json").write_text(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
