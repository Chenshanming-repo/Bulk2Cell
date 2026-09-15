"""Estimate empirical genomic TES boundaries from mapped Iso-Seq evidence.

Coordinates are zero-based half-open, and a TES is reported as an interbase
genomic boundary. Consensus records describe cluster endpoints; raw records
describe individual FLNC endpoints and remain distinct evidence types.
"""

from __future__ import annotations

import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping

from .models import Transcript


def _one_path(root: Path, patterns: tuple[str, ...]) -> Path | None:
    """Return the unique recursively discovered file matching ordered patterns."""
    for pattern in patterns:
        matches = sorted(path for path in root.rglob(pattern) if path.is_file())
        if len(matches) == 1:
            return matches[0]
        if len(matches) > 1:
            raise ValueError(f"ambiguous Iso-Seq files for {pattern}: {len(matches)} matches")
    return None


def discover_isoseq_paths(root: str | Path) -> dict[str, Path]:
    """Discover standard collapse assignment and mapped BAM files below a root."""
    base = Path(root)
    if not base.is_dir():
        raise ValueError(f"Iso-Seq root is not a directory: {base}")
    found: dict[str, Path] = {}
    group = _one_path(base, ("*collapsed.group.txt", "*collapse.group.txt"))
    read_stat = _one_path(base, ("*collapsed.read_stat.txt", "*collapse.read_stat.txt"))
    bam_candidates = sorted(path for path in base.rglob("*mapped.bam") if path.is_file())
    aligned = [path for path in bam_candidates if "04_align" in path.parts]
    if len(aligned) == 1:
        bam = aligned[0]
    elif len(bam_candidates) == 1:
        bam = bam_candidates[0]
    elif len(bam_candidates) > 1:
        raise ValueError(f"ambiguous consensus mapped BAM: {len(bam_candidates)} matches")
    else:
        bam = _one_path(base, ("*.bam",))
    if group:
        found["group"] = group
    if read_stat:
        found["read_stat"] = read_stat
    if bam:
        found["mapped_bam"] = bam
        found["consensus_bam"] = bam
    return found


def _record_assignment(
    assignments: dict[str, str], source_id: str, pbid: str
) -> None:
    """Store one assignment, rejecting an ID assigned to conflicting isoforms."""
    source_id, pbid = source_id.strip(), pbid.strip()
    if not source_id or not pbid:
        return
    previous = assignments.setdefault(source_id, pbid)
    if previous != pbid:
        raise ValueError(f"source ID {source_id} has conflicting PB.ID assignments")


def load_consensus_assignments(
    path: str | Path, selected_ids: Iterable[str] | None = None
) -> dict[str, str]:
    """Map collapsed consensus query IDs to retained Pigeon PB.IDs."""
    selected = set(selected_ids) if selected_ids is not None else None
    assignments: dict[str, str] = {}

    def rows():
        """Yield group member and PB.ID pairs without materializing the file."""
        handle = open(path)
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2:
                raise ValueError(f"{path}:{number}: expected PB.ID and group members")
            pbid = fields[0].strip()
            if pbid.lower() in {"pbid", "pb.id", "isoform"}:
                continue
            for source_id in fields[1].split(","):
                yield source_id, pbid
        handle.close()

    for source_id, pbid in rows():
        if selected is None or pbid in selected:
            _record_assignment(assignments, source_id, pbid)
    if selected is not None:
        # A second streaming pass catches an outside-PB conflict without retaining
        # millions of irrelevant member IDs in memory.
        for source_id, pbid in rows():
            if source_id.strip() in assignments and assignments[source_id.strip()] != pbid.strip():
                raise ValueError(f"source ID {source_id.strip()} has conflicting PB.ID assignments")
    return assignments


def load_raw_assignments(
    path: str | Path, selected_ids: Iterable[str] | None = None
) -> dict[str, str]:
    """Map raw FLNC IDs to retained Pigeon PB.IDs from read_stat."""
    selected = set(selected_ids) if selected_ids is not None else None
    assignments: dict[str, str] = {}

    def rows():
        """Yield validated raw ID and PB.ID pairs in a bounded-memory pass."""
        handle = open(path)
        reader = csv.DictReader(handle, delimiter="\t")
        names = {name.lower(): name for name in (reader.fieldnames or ())}
        id_col = names.get("id") or names.get("read_id")
        pb_col = names.get("pbid") or names.get("pb.id")
        length_col = names.get("length")
        if not id_col or not pb_col or not length_col:
            raise ValueError("read_stat requires id, length, and pbid columns")
        for number, row in enumerate(reader, 2):
            try:
                length = int(row[length_col])
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{path}:{number}: invalid read length") from exc
            if length <= 0:
                raise ValueError(f"{path}:{number}: read length must be positive")
            yield row[id_col], row[pb_col]
        handle.close()

    for source_id, pbid in rows():
        if selected is None or pbid in selected:
            _record_assignment(assignments, source_id, pbid)
    if selected is not None:
        for source_id, pbid in rows():
            if source_id.strip() in assignments and assignments[source_id.strip()] != pbid.strip():
                raise ValueError(f"source ID {source_id.strip()} has conflicting PB.ID assignments")
    return assignments


def _alignment_shape(read) -> tuple[tuple[tuple[int, int], ...], tuple[tuple[int, int], ...]]:
    """Return aligned reference blocks and explicit splice junctions from CIGAR."""
    position = read.reference_start
    blocks: list[tuple[int, int]] = []
    junctions: list[tuple[int, int]] = []
    for operation, length in read.cigartuples or ():
        if operation in (0, 7, 8):
            blocks.append((position, position + length))
            position += length
        elif operation == 3:
            junctions.append((position, position + length))
            position += length
        elif operation == 2:
            position += length
        elif operation not in (1, 4, 5, 6):
            raise ValueError(f"unsupported CIGAR operation {operation}")
    return tuple(blocks), tuple(junctions)


def _three_prime_softclip(read, strand: str) -> int:
    """Return soft- and hard-clipped bases at the transcript three-prime end."""
    cigar = read.cigartuples or ()
    terminal_ops = reversed(cigar) if strand == "+" else iter(cigar)
    clipped = 0
    for operation, length in terminal_ops:
        if operation not in (4, 5):
            break
        clipped += length
    return clipped


def _compatible(read, transcript: Transcript, tolerance: int, tes_window: int) -> tuple[bool, str, int]:
    """Validate locus, terminal-exon overlap, splice chain, and TES proximity."""
    if read.reference_name != transcript.chrom:
        return False, "wrong_locus", -1
    if bool(read.is_reverse) != (transcript.strand == "-"):
        return False, "wrong_strand", -1
    blocks, junctions = _alignment_shape(read)
    if not blocks:
        return False, "no_aligned_blocks", -1
    annotated = tuple(
        (transcript.exons[index][1], transcript.exons[index + 1][0])
        for index in range(len(transcript.exons) - 1)
    )
    # Internal junctions define identity; terminal exon endpoints may shift.
    if len(junctions) != len(annotated) or any(
        abs(observed[0] - expected[0]) > tolerance
        or abs(observed[1] - expected[1]) > tolerance
        for observed, expected in zip(junctions, annotated)
    ):
        return False, "incompatible_splice_chain", -1
    terminal_exon = transcript.exons[-1] if transcript.strand == "+" else transcript.exons[0]
    # CIGAR I/D and =/X boundaries split aligned blocks inside one exon.
    # Consider every aligned-base block in the terminal N-delimited exon;
    # neither a deletion span nor overlap from another exon supplies support.
    terminal_blocks = blocks
    if junctions:
        if transcript.strand == "+":
            terminal_blocks = tuple(block for block in blocks if block[0] >= junctions[-1][1])
        else:
            terminal_blocks = tuple(block for block in blocks if block[1] <= junctions[0][0])
    if not any(block[1] > terminal_exon[0] and block[0] < terminal_exon[1]
               for block in terminal_blocks):
        return False, "missing_terminal_exon", -1
    tes = read.reference_end if transcript.strand == "+" else read.reference_start
    if abs(tes - transcript.tes) > tes_window:
        return False, "tes_outside_window", tes
    return True, "", tes


def estimate_tes(
    bam_path: str | Path,
    transcripts: Iterable[Transcript],
    assignments: Mapping[str, str],
    *,
    mode: str,
    aliases: Mapping[str, str] | None = None,
    selected_ids: Iterable[str] | None = None,
    weight_mode: str = "consensus",
    min_mapq: int = 20,
    max_3prime_softclip: int = 20,
    junction_tolerance: int = 10,
    tes_window: int = 200,
) -> tuple[list[dict], list[dict], dict]:
    """Extract filtered endpoint records, per-isoform summaries, and QC.

    Raw mode counts each accepted read ID once. Consensus weighting gives every
    cluster record weight one by default. Support weighting instead uses the is
    tag; it cannot recover within-cluster endpoint variation.
    """
    import pysam

    if mode not in {"raw", "consensus"}:
        raise ValueError("mode must be raw or consensus")
    if weight_mode not in {"consensus", "support"}:
        raise ValueError("weight_mode must be consensus or support")
    if mode == "raw" and weight_mode != "consensus":
        raise ValueError("raw evidence must count each read once")
    if min(min_mapq, max_3prime_softclip, junction_tolerance, tes_window) < 0:
        raise ValueError("TES filter parameters must be nonnegative")
    models = {transcript.id: transcript for transcript in transcripts}
    alias_map = dict(aliases or {})
    selected = set(selected_ids) if selected_ids is not None else None
    evidence = "raw_read" if mode == "raw" else "consensus"
    qc = Counter()
    accepted: dict[str, tuple[str, int, float]] = {}
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(until_eof=True):
            qc["alignment_records"] += 1
            source_id = read.query_name
            assigned = assignments.get(source_id)
            if assigned is None:
                qc["unassigned_alignment"] += 1
                continue
            canonical = alias_map.get(assigned, assigned)
            if selected is not None and assigned not in selected and canonical not in selected:
                qc["unselected_isoform"] += 1
                continue
            transcript = models.get(canonical) or models.get(assigned)
            if transcript is None:
                qc["missing_transcript_model"] += 1
                continue
            if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail:
                qc["nonprimary"] += 1
                continue
            if read.mapping_quality < min_mapq:
                qc["low_mapq"] += 1
                continue
            if _three_prime_softclip(read, transcript.strand) > max_3prime_softclip:
                qc["excess_3prime_clip"] += 1
                continue
            compatible, reason, tes = _compatible(read, transcript, junction_tolerance, tes_window)
            if not compatible:
                qc[reason] += 1
                continue
            weight = 1.0
            if mode == "consensus" and weight_mode == "support":
                try:
                    if read.has_tag("is"):
                        weight = float(read.get_tag("is"))
                    else:
                        qc["missing_support_tag"] += 1
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"invalid is support tag for {source_id}") from exc
            if not math.isfinite(weight) or weight < 0:
                raise ValueError(f"invalid support weight for {source_id}")
            if weight == 0:
                qc["zero_weight"] += 1
                continue
            endpoint = (transcript.id, tes, weight)
            if source_id in accepted:
                if accepted[source_id] != endpoint:
                    raise ValueError(f"{source_id} has conflicting eligible alignments")
                qc["duplicate_read_id"] += 1
                continue
            accepted[source_id] = endpoint
            qc["accepted_alignments"] += 1
    weighted: dict[tuple[str, int], float] = defaultdict(float)
    observations: Counter[str] = Counter()
    for transcript_id, tes, weight in accepted.values():
        weighted[(transcript_id, tes)] += weight
        observations[transcript_id] += 1
    rows = [
        {"transcript_id": transcript_id, "tes": tes, "count": weight, "evidence": evidence}
        for (transcript_id, tes), weight in sorted(weighted.items())
    ]
    totals: dict[str, float] = defaultdict(float)
    for row in rows:
        totals[row["transcript_id"]] += row["count"]
    summary = [
        {"transcript_id": transcript_id, "observations": observations[transcript_id],
         "total_weight": totals[transcript_id], "evidence": evidence}
        for transcript_id in sorted(totals)
    ]
    result_qc = dict(qc)
    result_qc.update(
        evidence=evidence,
        weight_mode=weight_mode,
        parameters={"min_mapq": min_mapq, "max_3prime_softclip": max_3prime_softclip,
                    "junction_tolerance": junction_tolerance, "tes_window": tes_window},
        endpoint_rows=len(rows),
        isoforms=len(summary),
    )
    return rows, summary, result_qc


def write_tes_outputs(
    output_dir: str | Path,
    endpoints: Iterable[Mapping[str, object]],
    summary: Iterable[Mapping[str, object]],
    qc: Mapping[str, object],
) -> dict[str, Path]:
    """Write tes.tsv, summary.tsv, and completion metadata without overwrite."""
    output = Path(output_dir)
    paths = {"tes": output / "tes.tsv", "summary": output / "summary.tsv", "run": output / "run.json"}
    if any(path.exists() for path in paths.values()):
        raise FileExistsError(f"refusing to overwrite TES output in {output}")
    output.mkdir(parents=True, exist_ok=True)
    with paths["tes"].open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("transcript_id", "tes", "count", "evidence"), delimiter="\t")
        writer.writeheader()
        writer.writerows(endpoints)
    with paths["summary"].open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("transcript_id", "observations", "total_weight", "evidence"), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)
    payload = dict(qc)
    payload.update(status="complete", tes_path=str(paths["tes"]), summary_path=str(paths["summary"]))
    paths["run"].write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return paths
