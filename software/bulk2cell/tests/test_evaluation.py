"""Tests for sparse clustering and abundance evaluation utilities."""

from pathlib import Path

import numpy as np
import pytest
from scipy import sparse
from scipy.io import mmwrite

from bulk2cell.evaluation import (
    abundance_accuracy,
    compare_cell_matrices,
    evaluate_cells,
    load_matrix_market,
    read_cell_ranger_clusters,
)
import bulk2cell.evaluation as evaluation


def test_load_matrix_market_preserves_explicit_cell_feature_order(tmp_path: Path) -> None:
    """Load matrix market preserves explicit cell feature order."""
    matrix = sparse.csr_matrix([[1, 0, 2], [0, 3, 0]], dtype=float)
    path = tmp_path / "matrix.mtx"
    mmwrite(path, matrix)

    loaded, cells, features = load_matrix_market(
        path, ["cell-b", "cell-a"], ["tx-2", "tx-1", "tx-3"]
    )

    assert sparse.isspmatrix_csr(loaded)
    np.testing.assert_array_equal(loaded.toarray(), matrix.toarray())
    assert cells == ("cell-b", "cell-a")
    assert features == ("tx-2", "tx-1", "tx-3")


def test_load_matrix_market_rejects_shape_duplicate_ids_and_invalid_values(tmp_path: Path) -> None:
    """Load matrix market rejects shape duplicate ids and invalid values."""
    path = tmp_path / "matrix.mtx"
    mmwrite(path, sparse.csr_matrix([[1, 2], [3, 4]]))
    with pytest.raises(ValueError, match="shape"):
        load_matrix_market(path, ["a"], ["x", "y"])
    with pytest.raises(ValueError, match="unique"):
        load_matrix_market(path, ["a", "a"], ["x", "y"])

    bad = tmp_path / "bad.mtx"
    mmwrite(bad, sparse.csr_matrix([[1, -1], [0, np.inf]]))
    with pytest.raises(ValueError, match="nonnegative and finite"):
        load_matrix_market(bad, ["a", "b"], ["x", "y"])


def test_cell_ranger_csv_and_clustering_metrics_align_shuffled_barcodes(tmp_path: Path) -> None:
    """Cell ranger csv and clustering metrics align shuffled barcodes."""
    labels = tmp_path / "clusters.csv"
    labels.write_text("Barcode,Cluster\nc3,1\nc1,0\nc4,1\nc2,0\n")
    reference = read_cell_ranger_clusters(labels)
    matrix = sparse.csr_matrix(
        [[8, 1, 0], [0, 1, 8], [9, 0, 0], [0, 0, 9]], dtype=float
    )
    result = evaluate_cells(
        matrix,
        ["c1", "c4", "c2", "c3"],
        n_clusters=2,
        reference_labels=reference,
        n_components=2,
        random_state=7,
    )

    assert result.method == "TruncatedSVD"
    assert result.barcodes == ("c1", "c4", "c2", "c3")
    assert result.embedding.shape == (4, 2)
    assert result.adjusted_rand == pytest.approx(1.0)
    assert result.normalized_mutual_info == pytest.approx(1.0)
    assert np.isfinite(result.silhouette)
    np.testing.assert_array_equal(
        result.clusters,
        evaluate_cells(matrix, ["c1", "c4", "c2", "c3"], 2, random_state=7).clusters,
    )


def test_evaluate_cells_handles_sparse_zero_rows_and_clear_degenerate_errors() -> None:
    """Evaluate cells handles sparse zero rows and clear degenerate errors."""
    matrix = sparse.csr_matrix([[0, 0], [4, 0], [0, 4]], dtype=float)
    result = evaluate_cells(matrix, ["z", "a", "b"], 2, random_state=2)
    assert result.embedding.shape[0] == 3
    assert np.isfinite(result.embedding).all()

    with pytest.raises(ValueError, match="at least two cells"):
        evaluate_cells(sparse.csr_matrix([[1, 2]]), ["a"], 1)
    with pytest.raises(ValueError, match="n_clusters"):
        evaluate_cells(sparse.eye(3), ["a", "b", "c"], 3)
    with pytest.raises(ValueError, match="nonnegative and finite"):
        evaluate_cells(sparse.csr_matrix([[1, -2], [0, 1]]), ["a", "b"], 1)
    with pytest.raises(ValueError, match="distinct embedded rows"):
        evaluate_cells(sparse.csr_matrix(np.ones((3, 2))), ["a", "b", "c"], 2)


def test_silhouette_uses_bounded_distance_blocks(monkeypatch: pytest.MonkeyPatch) -> None:
    """Silhouette uses bounded distance blocks."""
    original_cdist = evaluation.cdist
    observed_rows: list[int] = []

    def recording_cdist(left: np.ndarray, right: np.ndarray) -> np.ndarray:
        """Recording cdist."""
        observed_rows.append(left.shape[0])
        return original_cdist(left, right)

    monkeypatch.setattr(evaluation, "cdist", recording_cdist)
    embedding = np.column_stack((np.arange(1200), np.zeros(1200)))
    score = evaluation._silhouette(embedding, np.repeat([0, 1], 600))
    assert np.isfinite(score)
    assert max(observed_rows) <= 512


def test_compare_original_scalpel_matrix_aligns_shuffled_barcodes() -> None:
    """Compare original scalpel matrix aligns shuffled barcodes."""
    current = sparse.csr_matrix([[8, 0], [0, 8], [7, 0], [0, 7]])
    original = sparse.csr_matrix([[0, 9], [9, 0], [0, 7], [7, 0]])
    comparison = compare_cell_matrices(
        current,
        ["a", "b", "c", "d"],
        original,
        ["d", "a", "b", "c"],
        n_clusters=2,
        random_state=4,
    )
    assert comparison.barcodes == ("a", "b", "c", "d")
    assert comparison.adjusted_rand == pytest.approx(1.0)
    assert comparison.normalized_mutual_info == pytest.approx(1.0)


def test_abundance_accuracy_aligns_ids_and_reports_expected_errors() -> None:
    """Abundance accuracy aligns ids and reports expected errors."""
    metrics = abundance_accuracy(
        truth=[1.0, 3.0, 2.0],
        estimate=[2.0, 1.0, 3.0],
        truth_ids=["a", "b", "c"],
        estimate_ids=["c", "a", "b"],
    )
    assert metrics.ids == ("a", "b", "c")
    assert metrics.mae == pytest.approx(0.0)
    assert metrics.rmse == pytest.approx(0.0)
    assert metrics.pearson == pytest.approx(1.0)
    assert metrics.spearman == pytest.approx(1.0)
    assert metrics.jensen_shannon == pytest.approx(0.0)
    assert metrics.compositional_distance == pytest.approx(0.0)


def test_abundance_accuracy_validates_alignment_overlap_and_degeneracy() -> None:
    """Abundance accuracy validates alignment overlap and degeneracy."""
    with pytest.raises(ValueError, match="same identifier set"):
        abundance_accuracy([1], [1], ["a"], ["b"])
    with pytest.raises(ValueError, match="supplied together"):
        abundance_accuracy([1], [1], ["tx"], ["tx"], training_read_ids=["r1"])
    with pytest.raises(ValueError, match="overlap"):
        abundance_accuracy(
            [1], [1], ["tx"], ["tx"],
            training_read_ids=["read-1"], validation_read_ids=["read-1"],
        )
    metrics = abundance_accuracy(
        [1], [1], ["tx"], ["tx"],
        training_read_ids=["train-read"], validation_read_ids=["heldout-read"],
    )
    assert metrics.mae == 0
    with pytest.raises(ValueError, match="nonnegative and finite"):
        abundance_accuracy([-1], [1], ["a"], ["a"])
    with pytest.raises(ValueError, match="positive total"):
        abundance_accuracy([0, 0], [0, 0], ["a", "b"], ["a", "b"])
