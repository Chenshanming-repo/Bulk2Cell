"""Evaluation utilities for bulk2cell cell matrices and abundances.

The cell workflow keeps barcode order explicit, library-size normalizes counts,
applies ``log1p``, embeds the sparse matrix with :class:`TruncatedSVD` (an
efficient PCA-like method that does not center sparse data), and clusters with
deterministically seeded K-means.  Cell Ranger clusters are reference labels,
not biological ground truth.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Hashable, Iterable, Mapping, Sequence

import numpy as np
from scipy import sparse
from scipy.io import mmread
from scipy.cluster.vq import kmeans2
from scipy.sparse.linalg import svds
from scipy.spatial.distance import cdist, jensenshannon
from scipy.stats import pearsonr, spearmanr


@dataclass(frozen=True)
class CellEvaluation:
    """Results from :func:`evaluate_cells`, in input barcode order."""

    barcodes: tuple[str, ...]
    clusters: np.ndarray
    embedding: np.ndarray
    method: str
    adjusted_rand: float | None
    normalized_mutual_info: float | None
    silhouette: float
    umap: np.ndarray | None = None


@dataclass(frozen=True)
class AbundanceMetrics:
    """ID-aligned scalar errors and distribution distances."""

    ids: tuple[Hashable, ...]
    mae: float
    rmse: float
    pearson: float
    spearman: float
    jensen_shannon: float
    compositional_distance: float


@dataclass(frozen=True)
class ClusterComparison:
    """Barcode-aligned agreement between clusterings of two cell matrices."""

    barcodes: tuple[str, ...]
    adjusted_rand: float
    normalized_mutual_info: float


def load_matrix_market(
    path: str | Path,
    cell_ids: Sequence[str],
    feature_ids: Sequence[str],
) -> tuple[sparse.csr_matrix, tuple[str, ...], tuple[str, ...]]:
    """Load a cell-by-feature Matrix Market file with explicit ordered IDs.

    Raises ``ValueError`` when dimensions disagree, IDs are duplicated, or
    values are negative/non-finite.  The matrix is returned in CSR format.
    """

    matrix = sparse.csr_matrix(mmread(path), dtype=float)
    cells = _unique_ids(cell_ids, "cell IDs")
    features = _unique_ids(feature_ids, "feature IDs")
    if matrix.shape != (len(cells), len(features)):
        raise ValueError(
            f"matrix shape {matrix.shape} does not match "
            f"{len(cells)} cell IDs and {len(features)} feature IDs"
        )
    _validate_matrix(matrix)
    return matrix, cells, features


def read_cell_ranger_clusters(path: str | Path) -> dict[str, str]:
    """Read a Cell Ranger CSV containing case-sensitive Barcode/Cluster columns."""

    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not {"Barcode", "Cluster"}.issubset(reader.fieldnames):
            raise ValueError("cluster CSV must contain Barcode and Cluster columns")
        result: dict[str, str] = {}
        for row in reader:
            barcode, cluster = row["Barcode"], row["Cluster"]
            if not barcode or barcode in result:
                raise ValueError("cluster CSV barcodes must be nonempty and unique")
            result[barcode] = cluster
    if not result:
        raise ValueError("cluster CSV contains no cells")
    return result


def evaluate_cells(
    matrix: sparse.spmatrix | np.ndarray,
    barcodes: Sequence[str],
    n_clusters: int,
    *,
    reference_labels: Mapping[str, object] | None = None,
    n_components: int = 20,
    random_state: int = 0,
    compute_umap: bool = False,
) -> CellEvaluation:
    """Normalize, embed, and cluster a cell-by-feature abundance matrix.

    Reference labels (for example Cell Ranger labels or labels derived from an
    original SCALPEL matrix) are aligned by barcode before ARI and NMI are
    computed.  All input barcodes must be present in the mapping.  ``umap`` is
    optional and raises an actionable ``ImportError`` when ``umap-learn`` is
    unavailable.  Zero-count cells remain valid zero rows.
    """

    x = sparse.csr_matrix(matrix, dtype=float)
    cells = _unique_ids(barcodes, "barcodes")
    if x.shape[0] != len(cells):
        raise ValueError("matrix row count must match barcode count")
    _validate_matrix(x)
    if x.shape[0] < 2 or x.shape[1] < 1:
        raise ValueError("evaluation requires at least two cells and one feature")
    if not isinstance(n_clusters, (int, np.integer)) or not 1 <= n_clusters < x.shape[0]:
        raise ValueError("n_clusters must be an integer from 1 to number of cells - 1")
    if n_components < 1:
        raise ValueError("n_components must be positive")

    totals = np.asarray(x.sum(axis=1)).ravel()
    scale = np.divide(1e4, totals, out=np.zeros_like(totals), where=totals > 0)
    normalized = sparse.diags(scale) @ x
    normalized.data = np.log1p(normalized.data)
    max_components = min(normalized.shape[0] - 1, normalized.shape[1])
    components = min(int(n_components), max_components)
    embedding = _svd_embedding(normalized, components, random_state)
    distinct_rows = np.unique(embedding, axis=0).shape[0]
    if n_clusters > distinct_rows:
        raise ValueError(
            f"n_clusters={n_clusters} exceeds the {distinct_rows} distinct embedded rows"
        )
    _, clusters = kmeans2(embedding, n_clusters, minit="++", seed=random_state)

    unique_clusters = np.unique(clusters)
    silhouette = (
        _silhouette(embedding, clusters)
        if 1 < len(unique_clusters) < len(cells)
        else float("nan")
    )
    ari: float | None = None
    nmi: float | None = None
    if reference_labels is not None:
        missing = [barcode for barcode in cells if barcode not in reference_labels]
        if missing:
            raise ValueError(f"reference labels missing {len(missing)} input barcodes")
        labels = [reference_labels[barcode] for barcode in cells]
        ari = _adjusted_rand(labels, clusters)
        nmi = _normalized_mutual_info(labels, clusters)

    umap_embedding = None
    if compute_umap:
        try:
            import umap
        except ImportError as exc:  # pragma: no cover - environment dependent
            raise ImportError("compute_umap=True requires the optional 'umap-learn' package") from exc
        umap_embedding = umap.UMAP(random_state=random_state).fit_transform(embedding)

    return CellEvaluation(cells, clusters, embedding, "TruncatedSVD", ari, nmi, silhouette, umap_embedding)


def abundance_accuracy(
    truth: Sequence[float],
    estimate: Sequence[float],
    truth_ids: Sequence[Hashable],
    estimate_ids: Sequence[Hashable],
    *,
    training_read_ids: Iterable[Hashable] | None = None,
    validation_read_ids: Iterable[Hashable] | None = None,
) -> AbundanceMetrics:
    """Compare truth/held-out abundances after exact identifier alignment.

    Training and validation read IDs may be supplied together to guard against
    circular evaluation; overlap between those read sets is rejected.  These
    IDs are independent of the transcript IDs being aligned. Jensen-Shannon is the
    square-root divergence returned by SciPy.  ``compositional_distance`` is
    total variation distance (half the L1 distance) between normalized
    compositions, which remains defined in the presence of zeros.  Pearson or
    Spearman correlation is ``nan`` when either aligned vector is constant.
    """

    truth_array = _abundance_vector(truth, "truth")
    estimate_array = _abundance_vector(estimate, "estimate")
    tids = _unique_hashable_ids(truth_ids, "truth IDs")
    eids = _unique_hashable_ids(estimate_ids, "estimate IDs")
    if len(tids) != len(truth_array) or len(eids) != len(estimate_array):
        raise ValueError("abundance vector lengths must match their identifier lengths")
    if set(tids) != set(eids):
        raise ValueError("truth and estimate must have the same identifier set")
    if (training_read_ids is None) != (validation_read_ids is None):
        raise ValueError("training_read_ids and validation_read_ids must be supplied together")
    if training_read_ids is not None and validation_read_ids is not None:
        training_reads = _iterable_id_set(training_read_ids, "training read IDs")
        validation_reads = _iterable_id_set(validation_read_ids, "validation read IDs")
        overlap = training_reads.intersection(validation_reads)
        if overlap:
            raise ValueError(f"training/validation read overlap would make evaluation circular: {len(overlap)} IDs")

    estimate_by_id = dict(zip(eids, estimate_array))
    aligned = np.asarray([estimate_by_id[item] for item in tids], dtype=float)
    if truth_array.sum() <= 0 or aligned.sum() <= 0:
        raise ValueError("truth and estimate must each have a positive total")
    delta = aligned - truth_array
    truth_composition = truth_array / truth_array.sum()
    estimate_composition = aligned / aligned.sum()
    pearson = _correlation(pearsonr, truth_array, aligned)
    spearman = _correlation(spearmanr, truth_array, aligned)
    return AbundanceMetrics(
        tids,
        float(np.mean(np.abs(delta))),
        float(np.sqrt(np.mean(delta**2))),
        pearson,
        spearman,
        float(jensenshannon(truth_composition, estimate_composition)),
        float(0.5 * np.abs(truth_composition - estimate_composition).sum()),
    )


def compare_cell_matrices(
    matrix: sparse.spmatrix | np.ndarray,
    barcodes: Sequence[str],
    original_matrix: sparse.spmatrix | np.ndarray,
    original_barcodes: Sequence[str],
    n_clusters: int,
    *,
    n_components: int = 20,
    random_state: int = 0,
) -> ClusterComparison:
    """Compare clusters from bulk2cell and an original SCALPEL matrix.

    Matrices may have different feature spaces.  Their cell identifiers must
    describe the same set; the original matrix is reordered to the current
    matrix's barcode order before deterministic, independent clustering.
    """

    cells = _unique_ids(barcodes, "barcodes")
    originals = _unique_ids(original_barcodes, "original barcodes")
    if set(cells) != set(originals):
        raise ValueError("current and original matrices must have the same barcode set")
    original = sparse.csr_matrix(original_matrix, dtype=float)
    if original.shape[0] != len(originals):
        raise ValueError("original matrix row count must match original barcode count")
    positions = {barcode: index for index, barcode in enumerate(originals)}
    aligned_original = original[[positions[barcode] for barcode in cells]]
    current_result = evaluate_cells(
        matrix, cells, n_clusters, n_components=n_components, random_state=random_state
    )
    original_result = evaluate_cells(
        aligned_original, cells, n_clusters, n_components=n_components, random_state=random_state
    )
    return ClusterComparison(
        cells,
        _adjusted_rand(current_result.clusters, original_result.clusters),
        _normalized_mutual_info(current_result.clusters, original_result.clusters),
    )


def _validate_matrix(matrix: sparse.csr_matrix) -> None:
    """Validate dimensions and stored sparse values without densifying."""

    if matrix.ndim != 2 or matrix.shape[0] == 0 or matrix.shape[1] == 0:
        raise ValueError("matrix must be a nonempty two-dimensional cell-by-feature matrix")
    if matrix.data.size and (not np.isfinite(matrix.data).all() or np.any(matrix.data < 0)):
        raise ValueError("matrix values must be nonnegative and finite")


def _svd_embedding(matrix: sparse.csr_matrix, components: int, seed: int) -> np.ndarray:
    """Return a deterministic PCA-like sparse TruncatedSVD embedding."""

    minimum_dimension = min(matrix.shape)
    if minimum_dimension == 1 or components >= minimum_dimension:
        _, _, right = np.linalg.svd(matrix.toarray(), full_matrices=False)
        return np.asarray(matrix @ right[:components].T)
    rng = np.random.default_rng(seed)
    _, singular, right = svds(matrix, k=components, v0=rng.standard_normal(minimum_dimension))
    order = np.argsort(singular)[::-1]
    return np.asarray(matrix @ right[order].T)


def _contingency(left: Sequence[object], right: Sequence[object]) -> np.ndarray:
    """Build a contingency table for two aligned categorical label vectors."""

    _, left_codes = np.unique(np.asarray(left, dtype=object), return_inverse=True)
    _, right_codes = np.unique(np.asarray(right, dtype=object), return_inverse=True)
    table = np.zeros((left_codes.max() + 1, right_codes.max() + 1), dtype=np.int64)
    np.add.at(table, (left_codes, right_codes), 1)
    return table


def _adjusted_rand(left: Sequence[object], right: Sequence[object]) -> float:
    """Compute adjusted Rand index from a contingency table."""

    table = _contingency(left, right)
    choose2 = lambda x: x * (x - 1) / 2
    observed = choose2(table).sum()
    rows = choose2(table.sum(axis=1)).sum()
    columns = choose2(table.sum(axis=0)).sum()
    total = choose2(table.sum())
    expected = rows * columns / total if total else 0.0
    denominator = (rows + columns) / 2 - expected
    return float((observed - expected) / denominator) if denominator else 1.0


def _normalized_mutual_info(left: Sequence[object], right: Sequence[object]) -> float:
    """Compute arithmetic-mean normalized mutual information."""

    table = _contingency(left, right).astype(float)
    probabilities = table / table.sum()
    row = probabilities.sum(axis=1)
    column = probabilities.sum(axis=0)
    rows, columns = np.nonzero(probabilities)
    mutual = np.sum(probabilities[rows, columns] * np.log(probabilities[rows, columns] / (row[rows] * column[columns])))
    h_left = -np.sum(row[row > 0] * np.log(row[row > 0]))
    h_right = -np.sum(column[column > 0] * np.log(column[column > 0]))
    denominator = (h_left + h_right) / 2
    return float(mutual / denominator) if denominator else 1.0


def _silhouette(embedding: np.ndarray, labels: np.ndarray) -> float:
    """Compute exact Euclidean silhouette using bounded row blocks."""

    values = np.zeros(len(labels), dtype=float)
    unique_labels = np.unique(labels)
    masks = {label: labels == label for label in unique_labels}
    for start in range(0, len(labels), 512):
        stop = min(start + 512, len(labels))
        distances = cdist(embedding[start:stop], embedding)
        for offset, index in enumerate(range(start, stop)):
            label = labels[index]
            same = masks[label]
            same_count = int(same.sum()) - 1
            if same_count == 0:
                continue
            within = distances[offset, same].sum() / same_count
            between = min(
                distances[offset, masks[other]].mean()
                for other in unique_labels
                if other != label
            )
            scale = max(within, between)
            values[index] = (between - within) / scale if scale else 0.0
    return float(values.mean())


def _unique_ids(values: Sequence[str], name: str) -> tuple[str, ...]:
    """Return validated nonempty, unique string identifiers."""

    result = tuple(values)
    if any(not isinstance(value, str) or not value for value in result):
        raise ValueError(f"{name} must be nonempty strings")
    if len(set(result)) != len(result):
        raise ValueError(f"{name} must be unique")
    return result


def _unique_hashable_ids(values: Sequence[Hashable], name: str) -> tuple[Hashable, ...]:
    """Return unique hashable identifiers while preserving order."""

    result = tuple(values)
    try:
        unique = set(result)
    except TypeError as exc:
        raise ValueError(f"{name} must be hashable") from exc
    if len(unique) != len(result):
        raise ValueError(f"{name} must be unique")
    return result


def _iterable_id_set(values: Iterable[Hashable], name: str) -> set[Hashable]:
    """Materialize a read-ID iterable as a set and validate hashability."""

    try:
        return set(values)
    except TypeError as exc:
        raise ValueError(f"{name} must be hashable") from exc


def _abundance_vector(values: Sequence[float], name: str) -> np.ndarray:
    """Convert and validate a one-dimensional abundance vector."""

    result = np.asarray(values, dtype=float)
    if result.ndim != 1 or result.size == 0:
        raise ValueError(f"{name} must be a nonempty one-dimensional vector")
    if not np.isfinite(result).all() or np.any(result < 0):
        raise ValueError(f"{name} values must be nonnegative and finite")
    return result


def _correlation(function: object, left: np.ndarray, right: np.ndarray) -> float:
    """Compute a correlation, returning nan for constant inputs."""

    if left.size < 2 or np.ptp(left) == 0 or np.ptp(right) == 0:
        return float("nan")
    return float(function(left, right).statistic)  # type: ignore[operator]
