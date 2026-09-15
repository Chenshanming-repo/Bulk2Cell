"""Synthetic tests for empirical Iso-Seq TES extraction."""

from pathlib import Path

import pytest

from bulk2cell.models import Transcript
from bulk2cell.tes import (
    discover_isoseq_paths,
    estimate_tes,
    load_consensus_assignments,
    load_raw_assignments,
    write_tes_outputs,
)


def _bam(tmp_path: Path, reads: list[dict]) -> Path:
    """Write and index coordinate-sorted synthetic alignments."""
    pysam = pytest.importorskip("pysam")
    raw = tmp_path / "raw.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 2000}]}
    with pysam.AlignmentFile(raw, "wb", header=header) as out:
        for row in reads:
            cigar = row["cigar"]
            read = pysam.AlignedSegment()
            read.query_name = row["name"]
            read.query_sequence = "A" * sum(n for op, n in cigar if op in (0, 1, 4, 7, 8))
            read.flag = row.get("flag", 0)
            read.reference_id = 0
            read.reference_start = row["start"]
            read.mapping_quality = row.get("mapq", 60)
            read.cigartuples = cigar
            read.query_qualities = pysam.qualitystring_to_array("I" * len(read.query_sequence))
            if "is" in row:
                read.set_tag("is", row["is"])
            out.write(read)
    result = tmp_path / "reads.bam"
    pysam.sort("-o", str(result), str(raw))
    pysam.index(str(result))
    return result


def test_assignment_parsers_select_pigeon_ids_and_detect_conflicts(tmp_path):
    """Consensus and raw tables map source IDs to one selected PB isoform."""
    group = tmp_path / "collapsed.group.txt"
    group.write_text("PB.1.1\tm640/a,transcript/2\nPB.2.1\tm640/b\n")
    assert load_consensus_assignments(group, {"PB.1.1"}) == {
        "m640/a": "PB.1.1",
        "transcript/2": "PB.1.1",
    }
    stat = tmp_path / "collapsed.read_stat.txt"
    stat.write_text("id\tlength\tpbid\nmovie/1/ccs\t800\tPB.1.1\nmovie/2/ccs\t700\tPB.2.1\n")
    assert load_raw_assignments(stat, {"PB.1.1"}) == {"movie/1/ccs": "PB.1.1"}
    stat.write_text("id\tlength\tpbid\na\t10\tPB.1.1\na\t10\tPB.2.1\n")
    with pytest.raises(ValueError, match="conflicting"):
        load_raw_assignments(stat)
    with pytest.raises(ValueError, match="conflicting"):
        load_raw_assignments(stat, {"PB.1.1"})
    stat.write_text(
        "id\tlength\tpbid\n"
        "outside-a\t10\tPB.9.1\nselected\t10\tPB.1.1\noutside-b\t10\tPB.8.1\n"
    )
    assert load_raw_assignments(stat, {"PB.1.1"}) == {"selected": "PB.1.1"}
    stat.write_text("id\tlength\tpbid\nselected\t10\tPB.9.1\nselected\t10\tPB.1.1\n")
    with pytest.raises(ValueError, match="conflicting"):
        load_raw_assignments(stat, {"PB.1.1"})


def test_consensus_endpoints_are_stranded_weighted_and_terminal_shift_tolerant(tmp_path):
    """Consensus support uses the is tag and accepts a shifted terminal boundary."""
    plus = Transcript("PB.1.1", "G1", "chr1", "+", ((100, 150), (200, 250)))
    minus = Transcript("PB.2.1", "G2", "chr1", "-", ((500, 550), (600, 650)))
    bam = _bam(tmp_path, [
        {"name": "c1", "start": 100, "cigar": [(0, 50), (3, 50), (0, 65)], "is": 4},
        {"name": "c3", "start": 100, "cigar": [(0, 50), (3, 50), (0, 65)]},
        {"name": "c2", "start": 485, "cigar": [(0, 65), (3, 50), (0, 50)], "flag": 16, "is": 3},
    ])
    rows, summary, qc = estimate_tes(
        bam, [plus, minus], {"c1": "PB.1.1", "c2": "PB.2.1", "c3": "PB.1.1"},
        mode="consensus", weight_mode="support",
    )
    assert [(r["transcript_id"], r["tes"], r["count"], r["evidence"]) for r in rows] == [
        ("PB.1.1", 265, 5.0, "consensus"),
        ("PB.2.1", 485, 3.0, "consensus"),
    ]
    assert {r["transcript_id"]: r["total_weight"] for r in summary} == {"PB.1.1": 5.0, "PB.2.1": 3.0}
    assert qc["accepted_alignments"] == 3 and qc["missing_support_tag"] == 1


def test_raw_reads_count_once_and_filter_primary_mapq_clip_strand_and_chain(tmp_path):
    """Raw FLNC evidence counts molecules once after all alignment quality checks."""
    tx = Transcript("PB.1.1", "G1", "chr1", "+", ((100, 150), (200, 250)))
    bam = _bam(tmp_path, [
        {"name": "good", "start": 100, "cigar": [(0, 50), (3, 50), (0, 50)]},
        {"name": "good", "start": 100, "cigar": [(0, 50), (3, 50), (0, 50)]},
        {"name": "low", "start": 100, "cigar": [(0, 50), (3, 50), (0, 50)], "mapq": 19},
        {"name": "clip", "start": 100, "cigar": [(0, 50), (3, 50), (0, 50), (4, 21)]},
        {"name": "hardclip", "start": 100, "cigar": [(0, 50), (3, 50), (0, 50), (4, 5), (5, 16)]},
        {"name": "reverse", "start": 100, "cigar": [(0, 50), (3, 50), (0, 50)], "flag": 16},
        {"name": "secondary", "start": 100, "cigar": [(0, 50), (3, 50), (0, 50)], "flag": 256},
        {"name": "chain", "start": 100, "cigar": [(0, 50), (3, 61), (0, 39)]},
    ])
    assignments = {name: "PB.1.1" for name in ("good", "low", "clip", "hardclip", "reverse", "secondary", "chain")}
    rows, _, qc = estimate_tes(bam, [tx], assignments, mode="raw")
    assert rows == [{"transcript_id": "PB.1.1", "tes": 250, "count": 1.0, "evidence": "raw_read"}]
    assert qc["duplicate_read_id"] == 1
    assert qc["low_mapq"] == qc["wrong_strand"] == 1
    assert qc["excess_3prime_clip"] == 2
    assert qc["nonprimary"] == qc["incompatible_splice_chain"] == 1


def test_conflicting_duplicate_primary_endpoints_are_rejected(tmp_path):
    """One raw ID cannot silently contribute an arbitrary one of two TES calls."""
    tx = Transcript("PB.1.1", "G1", "chr1", "+", ((100, 160),))
    bam = _bam(tmp_path, [
        {"name": "same", "start": 100, "cigar": [(0, 50)]},
        {"name": "same", "start": 100, "cigar": [(0, 60)]},
    ])
    with pytest.raises(ValueError, match="conflicting eligible alignments"):
        estimate_tes(bam, [tx], {"same": "PB.1.1"}, mode="raw")


def test_aliases_selection_zero_weight_and_unassigned_metadata(tmp_path):
    """Aliases canonicalize PB IDs, zero support is rejected, and exclusions are reported."""
    tx = Transcript("R1", "G", "chr1", "+", ((100, 150),))
    bam = _bam(tmp_path, [
        {"name": "zero", "start": 100, "cigar": [(0, 50)], "is": 0},
        {"name": "unknown", "start": 100, "cigar": [(0, 50)], "is": 2},
    ])
    rows, summary, qc = estimate_tes(
        bam, [tx], {"zero": "PB.1.1"}, mode="consensus", weight_mode="support",
        aliases={"PB.1.1": "R1"}, selected_ids={"R1"}
    )
    assert rows == [] and summary == []
    assert qc["zero_weight"] == 1 and qc["unassigned_alignment"] == 1
    assert qc["parameters"]["tes_window"] == 200


def test_discovery_and_output_tables_are_explicit_and_refuse_overwrite(tmp_path):
    """Discovery recognizes standard files and outputs preserve evidence provenance."""
    root = tmp_path / "isoseq"
    root.mkdir()
    for name in ("sample.collapsed.group.txt", "sample.collapsed.read_stat.txt", "sample.mapped.bam"):
        (root / name).write_text("")
    found = discover_isoseq_paths(root)
    assert found["group"].name.endswith("collapsed.group.txt")
    assert found["read_stat"].name.endswith("collapsed.read_stat.txt")
    align = root / "04_align"
    align.mkdir()
    (align / "consensus.mapped.bam").write_text("")
    (root / "raw.mapped.bam").write_text("")
    found = discover_isoseq_paths(root)
    assert found["mapped_bam"] == align / "consensus.mapped.bam"
    assert "raw_bam" not in found
    out = tmp_path / "output"
    paths = write_tes_outputs(
        out,
        [{"transcript_id": "PB.1", "tes": 42, "count": 2.0, "evidence": "consensus"}],
        [{"transcript_id": "PB.1", "observations": 1, "total_weight": 2.0, "evidence": "consensus"}],
        {"accepted_alignments": 1},
    )
    endpoints = paths["tes"]
    assert endpoints.read_text().splitlines()[0] == "transcript_id\ttes\tcount\tevidence"
    with pytest.raises(FileExistsError):
        write_tes_outputs(out, [], [], {})


@pytest.mark.parametrize('strand', ['+', '-'])
@pytest.mark.parametrize('middle', [[(1, 1)], [(2, 1)], [(8, 1)]])
def test_terminal_exon_overlap_survives_cigar_segments(tmp_path, strand, middle):
    """Indels and mismatch runs must not discard aligned terminal-exon overlap."""
    tx = Transcript('t', 'g', 'chr1', strand, ((100, 200),))
    cigar = [(7, 100), *middle, (7, 20)]
    if strand == '-':
        cigar = list(reversed(cigar))
    start = 100 if strand == '+' else 79
    bam = _bam(tmp_path, [{'name': 'r', 'start': start, 'cigar': cigar,
                           'flag': 16 if strand == '-' else 0}])
    rows, _, qc = estimate_tes(bam, [tx], {'r': 't'}, mode='raw')
    assert qc.get('accepted_alignments') == 1
    assert rows[0]['tes'] == (220 + (middle[0][0] != 1) if strand == '+' else 79)


@pytest.mark.parametrize('strand', ['+', '-'])
def test_deleted_bases_do_not_supply_terminal_exon_overlap(tmp_path, strand):
    """A deletion spanning an exon supplies no aligned-base overlap."""
    tx = Transcript('t', 'g', 'chr1', strand, ((100, 200),))
    bam = _bam(tmp_path, [{'name': 'r', 'start': 90,
                           'cigar': [(0, 10), (2, 100), (0, 10)],
                           'flag': 16 if strand == '-' else 0}])
    rows, _, qc = estimate_tes(bam, [tx], {'r': 't'}, mode='raw')
    assert rows == []
    assert qc['missing_terminal_exon'] == 1


@pytest.mark.parametrize('strand,exons,start,cigar', [
    ('+', ((100, 110), (115, 120)), 100, [(0, 25), (3, 5), (0, 10)]),
    ('-', ((100, 105), (110, 130)), 80, [(0, 10), (3, 5), (0, 35)]),
])
def test_overlap_across_splice_does_not_supply_terminal_exon(tmp_path, strand, exons, start, cigar):
    """Only the terminal N-delimited exon can supply terminal overlap."""
    tx = Transcript('t', 'g', 'chr1', strand, exons)
    bam = _bam(tmp_path, [{'name': 'r', 'start': start, 'cigar': cigar,
                           'flag': 16 if strand == '-' else 0}])
    rows, _, qc = estimate_tes(bam, [tx], {'r': 't'}, mode='raw', junction_tolerance=15)
    assert rows == []
    assert qc['missing_terminal_exon'] == 1
