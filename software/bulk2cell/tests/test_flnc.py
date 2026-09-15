"""Tests for selective FLNC FASTA extraction from PacBio BAM files."""

from __future__ import annotations

import gzip
import struct

import pytest

pysam = pytest.importorskip("pysam")

from bulk2cell.flnc import FLNCExtractionError, extract_flnc_fasta


def _write_bam(path, records):
    """Write unmapped query records and return their BGZF virtual offsets."""
    offsets = []
    with pysam.AlignmentFile(path, "wb", header={"HD": {"VN": "1.6"}}) as bam:
        for name, sequence in records:
            offsets.append(bam.tell())
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = sequence
            read.flag = 4
            read.query_qualities = pysam.qualitystring_to_array("I" * len(sequence))
            bam.write(read)
    return offsets


def _write_pbi(path, names, offsets, *, version=0x00040000):
    """Write the PBI header and basic-data section used by the extractor."""
    holes = [int(name.split("/")[1]) for name in names]
    count = len(names)
    with gzip.open(path, "wb") as pbi:
        pbi.write(struct.pack("<4sIHI18s", b"PBI\x01", version, 0, count, b"\0" * 18))
        pbi.write(struct.pack(f"<{count}i", *([0] * count)))  # rgId
        pbi.write(struct.pack(f"<{count}i", *([0] * count)))  # qStart
        pbi.write(struct.pack(f"<{count}i", *([0] * count)))  # qEnd
        pbi.write(struct.pack(f"<{count}i", *holes))
        pbi.write(struct.pack(f"<{count}f", *([1.0] * count)))
        pbi.write(bytes(count))
        pbi.write(struct.pack(f"<{count}q", *offsets))


def _fasta_records(path):
    """Return a small FASTA file as ordered name-to-sequence pairs."""
    lines = path.read_text().splitlines()
    return [(lines[index][1:], lines[index + 1]) for index in range(0, len(lines), 2)]


def test_pbi_extracts_only_exact_requested_names(tmp_path):
    """PBI candidates are filtered by hole and then by the complete query name."""
    bam = tmp_path / "reads.bam"
    records = [
        ("m1/10/ccs", "ACGT"),
        ("m2/10/ccs", "TTAA"),
        ("m84207_260409_155928_s1/11/ccs/0_4", "GGCC"),
    ]
    offsets = _write_bam(bam, records)
    pbi = tmp_path / "reads.bam.pbi"
    _write_pbi(pbi, [name for name, _ in records], offsets)
    output = tmp_path / "selected.fa"

    qc = extract_flnc_fasta(
        bam, {"m1/10/ccs", "m84207_260409_155928_s1/11/ccs/0_4"}, output
    )

    assert _fasta_records(output) == [
        ("m1/10/ccs", "ACGT"),
        ("m84207_260409_155928_s1/11/ccs/0_4", "GGCC"),
    ]
    assert qc["selected"] == 2
    assert qc["found"] == 2
    assert qc["missing"] == []
    assert qc["mode"] == "pbi"
    assert qc["pbi_version"] == 0x00040000
    assert qc["elapsed_seconds"] >= 0


def test_missing_read_does_not_publish_output_unless_allowed(tmp_path):
    """An incomplete selection is removed unless partial output was requested."""
    bam = tmp_path / "reads.bam"
    _write_bam(bam, [("movie/1/ccs", "AAAA")])
    output = tmp_path / "selected.fa"

    with pytest.raises(FLNCExtractionError, match="missing 1"):
        extract_flnc_fasta(bam, {"movie/1/ccs", "movie/2/ccs"}, output)
    assert not output.exists()

    qc = extract_flnc_fasta(
        bam, {"movie/1/ccs", "movie/2/ccs"}, output, allow_missing=True
    )
    assert qc["mode"] == "sequential"
    assert qc["missing"] == ["movie/2/ccs"]
    assert _fasta_records(output) == [("movie/1/ccs", "AAAA")]


def test_unsupported_name_and_pbi_version_fall_back_to_scan(tmp_path):
    """Names without PacBio holes and unsupported PBI versions use a safe scan."""
    bam = tmp_path / "reads.bam"
    offsets = _write_bam(bam, [("plain-name", "ACAC"), ("movie/8/ccs", "TGTG")])
    pbi = tmp_path / "custom.pbi"
    _write_pbi(pbi, ["movie/8/ccs", "movie/9/ccs"], offsets, version=0x00050000)

    name_qc = extract_flnc_fasta(bam, {"plain-name"}, tmp_path / "name.fa", pbi_path=pbi)
    version_qc = extract_flnc_fasta(bam, {"movie/8/ccs"}, tmp_path / "version.fa", pbi_path=pbi)

    assert name_qc["mode"] == "sequential"
    assert name_qc["pbi_version"] is None
    assert version_qc["mode"] == "sequential"
    assert version_qc["pbi_version"] == 0x00050000


def test_stale_pbi_offset_is_rejected(tmp_path):
    """An index row whose offset resolves to another query is not trusted."""
    bam = tmp_path / "reads.bam"
    names = ["movie/1/ccs", "movie/2/ccs"]
    offsets = _write_bam(bam, [(names[0], "AAAA"), (names[1], "CCCC")])
    pbi = tmp_path / "reads.bam.pbi"
    _write_pbi(pbi, names, list(reversed(offsets)))

    with pytest.raises(FLNCExtractionError, match="PBI offset mismatch"):
        extract_flnc_fasta(bam, {names[0]}, tmp_path / "selected.fa")


def test_duplicate_query_name_is_rejected(tmp_path):
    """Multiple BAM records for a requested exact name are treated as ambiguous."""
    bam = tmp_path / "reads.bam"
    names = ["movie/7/ccs", "movie/7/ccs"]
    offsets = _write_bam(bam, [(names[0], "AAAA"), (names[1], "CCCC")])
    _write_pbi(tmp_path / "reads.bam.pbi", names, offsets)

    with pytest.raises(FLNCExtractionError, match="duplicate query name"):
        extract_flnc_fasta(bam, {names[0]}, tmp_path / "selected.fa")


def test_refuses_to_overwrite_existing_output(tmp_path):
    """A completed or unrelated output file is never overwritten."""
    bam = tmp_path / "reads.bam"
    _write_bam(bam, [("movie/1/ccs", "AAAA")])
    output = tmp_path / "selected.fa"
    output.write_text("keep me\n")

    with pytest.raises(FileExistsError):
        extract_flnc_fasta(bam, {"movie/1/ccs"}, output)
    assert output.read_text() == "keep me\n"


def test_empty_selection_creates_empty_fasta_without_opening_bam(tmp_path):
    """An empty requested set completes without scanning or opening the input."""
    output = tmp_path / "empty.fa"

    qc = extract_flnc_fasta(tmp_path / "does-not-exist.bam", set(), output)

    assert output.read_bytes() == b""
    assert qc["selected"] == qc["found"] == 0
    assert qc["missing"] == []
