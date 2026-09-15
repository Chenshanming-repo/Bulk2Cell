#!/usr/bin/env python3
"""Export upstream benchmark tables as cells-by-feature sparse matrices."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path


def write_table(path, header, rows):
    """Write a standard headered TSV with deterministic row order."""
    with Path(path).open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def _stream_matrix(path, estimates, cell_index, feature_index, feature_key,
                   row_count, column_count, entry_count):
    """Write a Matrix Market coordinate file without retaining nonzeros."""
    with gzip.open(path, "wt") as handle:
        handle.write("%%MatrixMarket matrix coordinate real general\n")
        handle.write("% cells by features; SCALPEL relative abundance times "
                     "positive-probability molecules\n")
        handle.write(f"{row_count} {column_count} {entry_count}\n")
        with Path(estimates).open() as source:
            for row in csv.DictReader(source, delimiter="\t"):
                value = float(row["estimated_count"])
                if not math.isfinite(value) or value < 0:
                    raise ValueError("Nonfinite or negative upstream estimate")
                handle.write(f"{cell_index[row['bc']]+1} "
                             f"{feature_index[row[feature_key]]+1} {value:.17g}\n")


def export(directory, barcode_path):
    """Stream matrices while preserving requested cells and collapsed members."""
    directory, barcode_path = Path(directory), Path(barcode_path)
    completion_manifest = directory/"summary.json"
    if not completion_manifest.is_file():
        raise FileNotFoundError(f"SCALPEL run is incomplete: missing {completion_manifest}")
    published = [directory/name for name in (
        "isoform_counts.mtx.gz", "gene_counts.mtx.gz", "barcodes.tsv",
        "isoforms.tsv", "genes.tsv", "membership.tsv")]
    if any(path.exists() for path in published):
        raise FileExistsError(next(path for path in published if path.exists()))
    aliases = json.loads((directory/"aliases.json").read_text())
    cells = sorted(set(Path(barcode_path).read_text().splitlines()))
    models = {}
    for path in sorted(directory.glob("*/exons.tsv")):
        with path.open() as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                tx, gene = aliases[row["transcript_name"]], aliases[row["gene_name"]]
                members = (row["collapsed"].split("_") if row["collapsed"] != "none"
                           else [row["transcript_name"]])
                models[tx] = (gene, tuple(aliases[member] for member in members))
    transcripts = sorted(models)
    genes = sorted({value[0] for value in models.values()})
    cell_index = {value:index for index, value in enumerate(cells)}
    tx_index = {value:index for index, value in enumerate(transcripts)}
    gene_index = {value:index for index, value in enumerate(genes)}
    entry_count = 0
    total = 0.0
    group_count = 0
    previous_group = None
    estimates = directory/"estimates.tsv"
    with estimates.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            value = float(row["estimated_count"])
            if not math.isfinite(value) or value < 0:
                raise ValueError("Nonfinite or negative upstream estimate")
            if row["bc"] not in cell_index:
                raise ValueError(f"Estimate contains unrequested barcode: {row['bc']}")
            if row["transcript_name"] not in tx_index:
                raise ValueError(f"Estimate contains unknown transcript: {row['transcript_name']}")
            if row["gene_name"] not in gene_index:
                raise ValueError(f"Estimate contains unknown gene: {row['gene_name']}")
            entry_count += 1
            total += value
            group = (row["bc"], row["gene_name"])
            if group != previous_group:
                group_count += 1
                previous_group = group
    temporary = {path: path.with_name(path.name+".tmp") for path in published}
    try:
        _stream_matrix(temporary[directory/"isoform_counts.mtx.gz"], estimates,
                       cell_index, tx_index, "transcript_name",
                       len(cells), len(transcripts), entry_count)
        _stream_matrix(temporary[directory/"gene_counts.mtx.gz"], estimates,
                       cell_index, gene_index, "gene_name",
                       len(cells), len(genes), entry_count)
        write_table(temporary[directory/"barcodes.tsv"], ["barcode"],
                    ((value,) for value in cells))
        write_table(temporary[directory/"isoforms.tsv"],
                    ["transcript_id", "gene_id"],
                    ((value, models[value][0]) for value in transcripts))
        write_table(temporary[directory/"genes.tsv"], ["gene_id"],
                    ((value,) for value in genes))
        write_table(temporary[directory/"membership.tsv"],
                    ["representative_transcript_id", "gene_id", "member_transcript_id"],
                    ((value, models[value][0], member)
                     for value in transcripts for member in models[value][1]))
        for destination in published:
            temporary[destination].replace(destination)
    finally:
        for path in temporary.values():
            path.unlink(missing_ok=True)
    return {"cells": len(cells), "representative_isoforms": len(transcripts),
            "genes": len(genes), "estimated_count_sum": total,
            "estimated_cell_gene_groups": group_count, "matrix_entries": entry_count,
            "export_mode": "two-pass streaming Matrix Market"}


def main():
    """Provide an explicit export command for already completed upstream runs."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    parser.add_argument("barcodes", type=Path)
    args = parser.parse_args()
    summary = export(args.directory, args.barcodes)
    temporary = args.directory/"matrix_summary.json.tmp"
    temporary.write_text(json.dumps(summary, indent=2))
    temporary.replace(args.directory/"matrix_summary.json")


if __name__ == "__main__":
    main()
