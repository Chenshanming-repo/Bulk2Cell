"""Selective extraction of PacBio FLNC reads from BAM into FASTA.

PacBio ``.pbi`` files store BGZF virtual offsets for every BAM record.  This
module reads only the mandatory PBI basic-data section, identifies rows by ZMW
hole number, and verifies the complete query name after seeking into the BAM.
It therefore avoids a full pass over very large FLNC BAMs while remaining safe
when different movies contain the same hole number.
"""

from __future__ import annotations

import gzip
import os
from pathlib import Path
import struct
import tempfile
import time
from typing import BinaryIO, Iterable


class FLNCExtractionError(RuntimeError):
    """Report an incomplete or internally inconsistent FLNC extraction."""


_PBI_HEADER = struct.Struct("<4sIHI18s")
_PBI_MAGIC = b"PBI\x01"
_SUPPORTED_PBI_MAJORS = {3, 4}


def _pbi_major(version: int) -> int:
    """Return the major component used by known PBI integer encodings."""
    high_major = (version >> 24) & 0xFF
    return high_major or ((version >> 16) & 0xFF)


def _read_exact(stream: BinaryIO, size: int, label: str) -> bytes:
    """Read exactly ``size`` bytes or report a truncated PBI section."""
    value = stream.read(size)
    if len(value) != size:
        raise FLNCExtractionError(f"truncated PBI {label}")
    return value


def _discard_exact(stream: BinaryIO, size: int, label: str) -> None:
    """Consume a PBI section in bounded chunks and detect truncation."""
    remaining = size
    while remaining:
        chunk_size = min(1_048_576, remaining)
        _read_exact(stream, chunk_size, label)
        remaining -= chunk_size


def _open_pbi(path: Path) -> BinaryIO:
    """Open a PBI file whether it is BGZF/gzip-compressed or plain binary."""
    raw = path.open("rb")
    signature = raw.read(2)
    raw.seek(0)
    if signature == b"\x1f\x8b":
        raw.close()
        return gzip.open(path, "rb")
    return raw


def _pacbio_hole(query_name: str) -> int | None:
    """Parse a hole from PacBio consensus or subread-style query names.

    Accepted forms are ``movie/hole/ccs``, ``movie/hole/start_end``, and the
    segmented-consensus form ``movie/hole/ccs/start_end`` used by HBA8 FLNC.
    """
    fields = query_name.split("/")
    if len(fields) not in {3, 4} or not fields[0] or not fields[1].isdigit():
        return None
    interval_field: str | None
    if len(fields) == 4:
        if fields[2] != "ccs":
            return None
        interval_field = fields[3]
    else:
        interval_field = None if fields[2] == "ccs" else fields[2]
    if interval_field is not None:
        interval = interval_field.split("_")
        if len(interval) != 2 or not all(value.isdigit() for value in interval):
            return None
    return int(fields[1])


def _read_pbi_offsets(path: Path, selected_holes: set[int]) -> tuple[int, list[tuple[int, int]]]:
    """Return PBI version and ``(hole, virtual_offset)`` candidate rows.

    Only matching row indexes and offsets are retained.  The full hole and
    offset columns are consumed in chunks to keep memory independent of the
    total number of BAM records.
    """
    with _open_pbi(path) as stream:
        magic, version, _flags, count, _reserved = _PBI_HEADER.unpack(
            _read_exact(stream, _PBI_HEADER.size, "header")
        )
        if magic != _PBI_MAGIC:
            raise FLNCExtractionError(f"invalid PBI magic in {path}")
        if _pbi_major(version) not in _SUPPORTED_PBI_MAJORS:
            return version, []

        # Skip rgId, qStart, and qEnd.  holeNumber is the fourth basic column.
        _discard_exact(stream, count * 12, "basic prefix")
        candidate_rows: dict[int, int] = {}
        for row_start in range(0, count, 65_536):
            chunk_count = min(65_536, count - row_start)
            chunk = _read_exact(stream, chunk_count * 4, "holeNumber column")
            for local_row, (hole,) in enumerate(struct.iter_unpack("<i", chunk)):
                if hole in selected_holes:
                    candidate_rows[row_start + local_row] = hole

        # readQual and ctxtFlag precede fileOffset in the basic data section.
        _discard_exact(stream, count * 5, "basic suffix")
        candidates: list[tuple[int, int]] = []
        for row_start in range(0, count, 32_768):
            chunk_count = min(32_768, count - row_start)
            chunk = _read_exact(stream, chunk_count * 8, "fileOffset column")
            for local_row, (offset,) in enumerate(struct.iter_unpack("<q", chunk)):
                row = row_start + local_row
                if row in candidate_rows:
                    candidates.append((candidate_rows[row], offset))
    return version, candidates


def _forward_sequence(read: object) -> str:
    """Return a BAM record's query sequence in its original forward orientation."""
    sequence = read.get_forward_sequence()
    if sequence is None:
        raise FLNCExtractionError(f"query {read.query_name!r} has no sequence")
    return sequence


def _extract_sequential(bam: object, selected: set[str], output: BinaryIO) -> set[str]:
    """Scan a BAM once, writing exact requested names and rejecting duplicates."""
    found: set[str] = set()
    for read in bam.fetch(until_eof=True):
        name = read.query_name
        if name not in selected:
            continue
        if name in found:
            raise FLNCExtractionError(f"duplicate query name in BAM: {name}")
        output.write(f">{name}\n{_forward_sequence(read)}\n")
        found.add(name)
    return found


def _extract_indexed(
    bam: object,
    selected: set[str],
    candidates: list[tuple[int, int]],
    output: BinaryIO,
) -> set[str]:
    """Seek to PBI candidates, verify their identities, and write selected reads."""
    found: set[str] = set()
    for indexed_hole, offset in candidates:
        try:
            bam.seek(offset)
            read = next(bam)
        except (OSError, StopIteration, ValueError) as error:
            raise FLNCExtractionError(f"invalid PBI virtual offset {offset}") from error
        actual_hole = _pacbio_hole(read.query_name)
        if actual_hole != indexed_hole:
            raise FLNCExtractionError(
                f"PBI offset mismatch: indexed hole {indexed_hole}, "
                f"found query {read.query_name!r}"
            )
        # Same hole numbers can occur in multiple movies, so only exact names
        # are emitted.  A nonmatching movie is a legitimate candidate collision.
        if read.query_name not in selected:
            continue
        if read.query_name in found:
            raise FLNCExtractionError(f"duplicate query name in BAM: {read.query_name}")
        output.write(f">{read.query_name}\n{_forward_sequence(read)}\n")
        found.add(read.query_name)
    return found


def extract_flnc_fasta(
    bam_path: str | os.PathLike[str],
    read_ids: Iterable[str],
    output_path: str | os.PathLike[str],
    pbi_path: str | os.PathLike[str] | None = None,
    *,
    allow_missing: bool = False,
) -> dict[str, object]:
    """Extract selected FLNC BAM queries as forward-orientation FASTA records.

    Existing outputs are refused.  When a supported adjacent ``.pbi`` exists
    and every requested name has a standard PacBio form, BAM records are read
    by virtual offset.  A missing PBI, an unsupported PBI version, or a
    nonstandard requested name triggers one sequential BAM scan.  Malformed or
    stale supported indexes raise :class:`FLNCExtractionError` rather than
    silently scanning the BAM.

    Output is first written beside the destination and atomically published
    only after validation.  By default any missing requested ID aborts and
    removes the temporary file; ``allow_missing=True`` explicitly permits a
    partial FASTA.  Returned QC contains selected/found counts, sorted missing
    IDs, access mode, raw PBI version, paths, and elapsed wall-clock seconds.
    """
    import pysam

    started = time.monotonic()
    bam_file = Path(bam_path)
    destination = Path(output_path)
    if destination.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {destination}")
    selected = {str(read_id) for read_id in read_ids}
    index_file = Path(pbi_path) if pbi_path is not None else Path(f"{bam_file}.pbi")
    holes = {_pacbio_hole(name) for name in selected}
    use_index = bool(selected) and None not in holes and index_file.is_file()
    pbi_version: int | None = None
    candidates: list[tuple[int, int]] = []
    if use_index:
        pbi_version, candidates = _read_pbi_offsets(index_file, {int(hole) for hole in holes})
        if _pbi_major(pbi_version) not in _SUPPORTED_PBI_MAJORS:
            use_index = False

    destination.parent.mkdir(parents=True, exist_ok=True)
    temp_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="ascii", newline="\n", dir=destination.parent,
            prefix=f".{destination.name}.", suffix=".tmp", delete=False,
        ) as output:
            temp_name = output.name
            if not selected:
                found = set()
                mode = "sequential"
            else:
                # Raw FLNC BAMs may contain only unmapped queries and no SQ header.
                with pysam.AlignmentFile(str(bam_file), "rb", check_sq=False) as bam:
                    if use_index:
                        found = _extract_indexed(bam, selected, candidates, output)
                        mode = "pbi"
                    else:
                        found = _extract_sequential(bam, selected, output)
                        mode = "sequential"
        missing = sorted(selected - found)
        if missing and not allow_missing:
            raise FLNCExtractionError(
                f"missing {len(missing)} of {len(selected)} requested FLNC reads"
            )
        # Hard-link publication is atomic and fails if a destination appeared
        # since the initial existence check, preserving the no-overwrite rule.
        try:
            os.link(temp_name, destination)
        except FileExistsError as error:
            raise FileExistsError(
                f"refusing to overwrite existing output: {destination}"
            ) from error
        os.unlink(temp_name)
        temp_name = None
    finally:
        if temp_name is not None:
            try:
                os.unlink(temp_name)
            except FileNotFoundError:
                pass

    return {
        "selected": len(selected),
        "found": len(found),
        "missing": missing,
        "mode": mode,
        "pbi_version": pbi_version,
        "bam_path": str(bam_file),
        "pbi_path": str(index_file) if index_file.is_file() else None,
        "output_path": str(destination),
        "elapsed_seconds": time.monotonic() - started,
    }
