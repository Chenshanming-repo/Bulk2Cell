"""Tests for the matched multi-method HBA8 benchmark evaluator."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
from scipy.io import mmwrite


def _module():
    """Load the standalone benchmark script as a module."""

    path = Path(__file__).parents[1] / "benchmark/evaluate_hba8.py"
    spec = importlib.util.spec_from_file_location("evaluate_hba8", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _write_matrix(path: Path, values: list[list[float]]) -> None:
    """Write a compressed Matrix Market fixture."""

    with gzip.open(path, "wb") as handle:
        mmwrite(handle, sparse.csr_matrix(values))


def _write_tsv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    """Write a small benchmark metadata fixture."""

    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def _run_fixture(root: Path, kind: str, matrix: list[list[float]]) -> Path:
    """Create one completed bulk2cell or SCALPEL-style run."""

    root.mkdir()
    _write_matrix(root / "isoform_counts.mtx.gz", matrix)
    _write_tsv(root / "barcodes.tsv", ["barcode"], [["c3"], ["c1"], ["c4"], ["c2"]])
    _write_tsv(
        root / "isoforms.tsv",
        ["transcript_id", "gene_id"],
        [["tx1", "g1"], ["tx2", "g1"], ["tx3", "g2"]],
    )
    if kind == "bulk2cell":
        (root / "run.json").write_text(json.dumps({
            "status": "complete", "elapsed_seconds": 12.5, "peak_rss_mib": 42.0,
            "cells": 4, "isoforms": 3, "genes": 2,
            "qc": {"assigned_molecules": 15, "incompatible_molecules": 5,
                   "em_nonconverged_cells": 2},
        }))
        _write_tsv(
            root / "gene_qc.tsv",
            ["gene_id", "isoforms", "groups", "input_molecules", "assigned_molecules",
             "incompatible_molecules", "em_nonconverged_cells", "max_em_iterations"],
            [["g1", 2, 2, 10, 8, 2, 1, 500], ["g2", 1, 1, 10, 7, 3, 1, 500]],
        )
    else:
        (root / "summary.json").write_text(json.dumps({
            "wall_seconds": 20.0, "peak_child_rss_kib": 10240,
            "convergence": "Upstream MAX_IT=30, tolerance=0.01; iteration counts not exposed",
            "count_scaling": "relative abundance times supported molecules",
        }))
        (root / "matrix_summary.json").write_text(json.dumps({"cells": 4, "genes": 2}))
    return root


def test_load_run_requires_a_true_completion_manifest(tmp_path: Path) -> None:
    """Directories without a successful method-specific manifest are rejected."""

    benchmark = _module()
    run = _run_fixture(tmp_path / "run", "bulk2cell", [[1, 0, 0]] * 4)
    manifest = json.loads((run / "run.json").read_text())
    manifest["status"] = "incomplete"
    (run / "run.json").write_text(json.dumps(manifest))

    with pytest.raises(ValueError, match="completed bulk2cell run"):
        benchmark.load_run("method", run, "bulk2cell")


def test_run_summary_aligns_barcodes_and_aggregates_gene_totals(tmp_path: Path) -> None:
    """Callability and gene totals are derived from explicit matrix metadata."""

    benchmark = _module()
    run = _run_fixture(
        tmp_path / "run", "bulk2cell",
        [[3, 1, 0], [0, 0, 0], [0, 0, 4], [2, 0, 1]],
    )
    loaded = benchmark.load_run("bulk", run, "bulk2cell")
    summary, gene_rows = benchmark.summarize_run(loaded)

    assert summary["whitelist_cells"] == 4
    assert summary["active_cells"] == 3
    assert summary["callable_cell_gene_pairs"] == 4
    assert summary["total_estimated_counts"] == pytest.approx(11.0)
    assert summary["runtime_seconds"] == pytest.approx(12.5)
    assert summary["peak_rss_mib"] == pytest.approx(42.0)
    assert [(row[1], row[2]) for row in gene_rows] == [("g1", 6.0), ("g2", 5.0)]


def test_shared_active_cells_are_intersected_by_barcode(tmp_path: Path) -> None:
    """Shared comparisons use the intersection of nonzero cells, independent of order."""

    benchmark = _module()
    left = benchmark.load_run(
        "left", _run_fixture(tmp_path / "left", "bulk2cell",
                             [[1, 0, 0], [1, 0, 0], [0, 0, 0], [1, 0, 0]]),
        "bulk2cell",
    )
    right = benchmark.load_run(
        "right", _run_fixture(tmp_path / "right", "scalpel",
                              [[1, 0, 0], [0, 0, 0], [1, 0, 0], [1, 0, 0]]),
        "scalpel",
    )

    assert benchmark.shared_active_barcodes([left, right]) == ("c2", "c3")


def test_evaluate_shared_writes_metrics_coordinates_and_standalone_plot(tmp_path: Path) -> None:
    """Shared-cell evaluation emits aligned tables, strict JSON, and a labeled PNG."""

    benchmark = _module()
    matrix = [[9, 0, 0], [0, 0, 9], [8, 0, 0], [0, 0, 8]]
    run = benchmark.load_run(
        "bulk", _run_fixture(tmp_path / "run", "bulk2cell", matrix), "bulk2cell"
    )
    labels = {"c1": "A", "c2": "A", "c3": "B", "c4": "B"}
    out = tmp_path / "evaluation"
    result = benchmark.evaluate_shared(
        run, ("c1", "c2", "c3", "c4"), labels, out,
        n_clusters=2, seed=3, compute_umap=False,
    )

    assert result["shared_active_cells"] == 4
    assert result["ari_cell_ranger_reference"] == pytest.approx(1.0)
    assert (out / "bulk_shared_embedding.tsv").exists()
    assert (out / "bulk_shared_clusters.tsv").exists()
    assert (out / "bulk_shared_svd.png").stat().st_size > 0
    assert "NaN" not in json.dumps(result, allow_nan=False)


def test_full_sample_scope_reports_actual_gene_universe(tmp_path: Path) -> None:
    """Full-sample reports derive gene counts and do not retain the old17-gene claim."""
    benchmark=_module()
    path=_run_fixture(tmp_path/'run','bulk2cell',[[9,0,0],[0,0,9],[8,0,0],[0,0,8]])
    labels=tmp_path/'clusters.csv';labels.write_text('Barcode,Cluster\nc1,A\nc2,A\nc3,B\nc4,B\n')
    out=tmp_path/'evaluation'
    assert benchmark.main(['--run',f'full:bulk2cell:{path}','--clusters',str(labels),
        '--out',str(out),'--scope','full-sample','--n-clusters','2','--no-umap'])==0
    report=json.loads((out/'evaluation.json').read_text())
    assert report['scope']=='full-sample'
    assert report['genes_union']==report['genes_shared']==2
    assert report['panel_genes']==2
    assert not any('17-gene' in value or 'targeted' in value for value in report['limitations'])


def test_shared_runner_resource_scope_is_preserved(tmp_path: Path) -> None:
    """Joint variant timings must not appear as independent method resources."""
    benchmark = _module()
    path = _run_fixture(tmp_path / 'run', 'bulk2cell', [[1, 0, 0]] * 4)
    manifest = json.loads((path / 'run.json').read_text())
    manifest.update(runtime_scope='shared_three_variant_runner',
                    peak_rss_scope='parent_only', external_training_elapsed_seconds=17.0,
                    fit_metrics={'variant_worker_seconds': {'empirical_tes': 8.0},
                                 'peak_rss_mib': 91.0})
    (path / 'run.json').write_text(json.dumps(manifest))
    summary, _ = benchmark.summarize_run(benchmark.load_run('bulk', path, 'bulk2cell'))
    assert summary['runtime_scope'] == 'shared_three_variant_runner'
    assert summary['peak_rss_scope'] == 'parent_only'
    assert summary['external_training_elapsed_seconds'] == 17.0
    assert summary['fit_metrics']['peak_rss_mib'] == 91.0
