"""Behavioral tests for bulk2cell input adapters."""

import json
from pathlib import Path

import pytest

from bulk2cell.adapters import (
    extract_molecules,
    load_barcodes,
    load_regionquant_catalog,
    load_salmon_quant,
    load_tes_tsv,
    load_transcript_catalog,
)


def test_reference_pigeon_union_deduplicates_structure_and_keeps_novel_genes(tmp_path):
    """Reference aliases aggregate LR support while a novel LR gene remains present."""
    gtf = tmp_path / "ref.gtf"
    gtf.write_text(
        'chr1\tr\texon\t11\t20\t.\t+\t.\tgene_id "ENSG1.4"; gene_name "A"; transcript_id "R1";\n'
        'chr1\tr\texon\t31\t40\t.\t+\t.\tgene_id "ENSG1.4"; gene_name "A"; transcript_id "R1";\n'
        'chr1\tr\texon\t11\t20\t.\t+\t.\tgene_id "ENSG1.4"; gene_name "A"; transcript_id "R_ALIAS";\n'
        'chr1\tr\texon\t31\t40\t.\t+\t.\tgene_id "ENSG1.4"; gene_name "A"; transcript_id "R_ALIAS";\n'
    )
    gff = tmp_path / "pigeon.gff"
    gff.write_text(
        "chr1\tp\texon\t11\t20\t.\t+\t.\ttranscript_id \"PB.1\";\n"
        "chr1\tp\texon\t31\t40\t.\t+\t.\ttranscript_id \"PB.1\";\n"
        "chr2\tp\texon\t101\t130\t.\t-\t.\tParent=PB.2\n"
    )
    cls = tmp_path / "class.tsv"
    cls.write_text(
        "isoform\tchrom\tstrand\tassociated_gene\tstructural_category\tfl_assoc\n"
        "PB.1\tchr1\t+\tENSG1.9\tfull-splice_match\t7\n"
        "PB.2\tchr2\t-\tNOVELG\tnovel_not_in_catalog\t3\n"
    )
    tx, qc = load_transcript_catalog(gtf, gff, cls)
    assert [(x.id, x.gene_id, x.exons, x.lr_count) for x in tx] == [
        ("R1", "ENSG1.4", ((10, 20), (30, 40)), 7.0),
        ("PB.2", "NOVELG", ((100, 130),), 3.0),
    ]
    assert {k for k, v in qc["aliases"].items() if v == "R1"} == {"R1", "R_ALIAS", "PB.1"}
    assert qc["deduplicated_transcripts"] == 2


def test_regionquant_and_abundance_adapters_aggregate_aliases(tmp_path):
    """Existing catalog coordinates and quantification aliases use the same structure model."""
    catalog = tmp_path / "catalog.json"
    catalog.write_text(json.dumps({"genes": {"ENSG1": {"gene_id": "ENSG1", "gene_name": "A", "chrom": "chr1", "strand": "+", "ref": [{"id": "R1", "exons": [[0, 10]]}], "lr": [{"id": "PB1", "exons": [[0, 10]], "fl": 4}]}}}))
    tx, qc = load_regionquant_catalog(catalog)
    assert len(tx) == 1 and tx[0].lr_count == 4
    quant = tmp_path / "quant.sf"
    quant.write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nR1\t10\t8\t2\t5\nPB1\t10\t8\t3\t7\n")
    values, qqc = load_salmon_quant(quant, tx, qc["aliases"])
    assert values == {tx[0].id: 5.0}
    assert qqc["matched_rows"] == 2


def test_tes_and_barcode_files_validate_and_preserve_identifiers(tmp_path):
    """TES observations are numeric and full GEM barcode suffixes are retained."""
    tes = tmp_path / "tes.tsv"
    tes.write_text("transcript_id\ttes\nR1\t40\nR1\t42\n")
    values, qc = load_tes_tsv(tes)
    assert values == {"R1": [(40, 1.0), (42, 1.0)]} and qc["observations"] == 2
    barcodes = tmp_path / "barcodes.tsv"
    barcodes.write_text("AAAC-1\nAAAG-1\nAAAC-1\n")
    assert load_barcodes(barcodes) == {"AAAC-1", "AAAG-1"}


def test_bam_extraction_unions_conflicting_read_evidence_and_counts_rejections(tmp_path):
    """A molecule retains the union of all eligible read blocks/junctions and reports GX exclusions."""
    pysam = pytest.importorskip("pysam")
    from bulk2cell.models import Transcript

    raw = tmp_path / "reads.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(raw, "wb", header=header) as out:
        def emit(name, start, cigar, tags, reverse=False):
            """Emit."""
            read = pysam.AlignedSegment()
            read.query_name = name
            read.query_sequence = "A" * sum(n for op, n in cigar if op in (0, 1, 4, 7, 8))
            read.flag = 16 if reverse else 0
            read.reference_id = 0; read.reference_start = start; read.mapping_quality = 255
            read.cigartuples = cigar; read.query_qualities = pysam.qualitystring_to_array("I" * len(read.query_sequence))
            read.set_tags(tags); out.write(read)
        common = [("CB", "AAAC-1"), ("UB", "U1"), ("GX", "ENSG1.8")]
        emit("a", 100, [(0, 10), (3, 10), (0, 10)], common)
        emit("b", 105, [(0, 10), (3, 20), (0, 5)], common)
        emit("missing", 100, [(0, 5)], [("CB", "AAAC-1"), ("UB", "U2")])
        emit("multi", 100, [(0, 5)], [("CB", "AAAC-1"), ("UB", "U3"), ("GX", "ENSG1;ENSG2")])
        emit("wrongstrand", 100, [(0, 5)], [("CB", "AAAC-1"), ("UB", "U4"), ("GX", "ENSG1")], True)
    bam = tmp_path / "reads.sorted.bam"
    pysam.sort("-o", str(bam), str(raw)); pysam.index(str(bam))
    transcript = Transcript("R1", "ENSG1.4", "chr1", "+", ((90, 200),))
    molecules, qc = extract_molecules(bam, [transcript], {"ENSG1.4"}, {"AAAC-1"})
    assert len(molecules) == 1
    assert molecules[0].gene_id == "ENSG1.4"
    assert molecules[0].blocks == ((100, 115), (120, 130), (135, 140))
    assert molecules[0].junctions == ((110, 120), (115, 135))
    assert qc["eligible_reads"] == 2
    assert qc["missing_gx"] == qc["multigene_gx"] == qc["wrong_strand"] == 1



def test_salmon_prefers_tpm_and_rejects_duplicate_or_nonfinite_rows(tmp_path):
    """TPM is the abundance measure; fallback and malformed input are explicit."""
    from bulk2cell.models import Transcript
    tx = [Transcript("R1", "G", "chr1", "+", ((0, 10),))]
    quant = tmp_path / "quant.sf"
    quant.write_text("Name\tTPM\tNumReads\nR1\t2.5\t99\n")
    values, qc = load_salmon_quant(quant, tx)
    assert values == {"R1": 2.5} and qc["measure"] == "TPM"
    quant.write_text("Name\tNumReads\nR1\t4\n")
    values, qc = load_salmon_quant(quant, tx)
    assert values == {"R1": 4.0} and qc["measure"] == "NumReads_fallback"
    quant.write_text("Name\tTPM\tNumReads\nR1\t1\t1\nR1\t2\t2\n")
    with pytest.raises(ValueError, match="duplicate"):
        load_salmon_quant(quant, tx)
    quant.write_text("Name\tTPM\tNumReads\nR1\tnan\t1\n")
    with pytest.raises(ValueError, match="finite"):
        load_salmon_quant(quant, tx)


def test_tes_rejects_fractional_nonfinite_and_zero_total_weights(tmp_path):
    """TES boundaries are integers and each transcript needs positive finite mass."""
    tes = tmp_path / "tes.tsv"
    tes.write_text("transcript_id\ttes\tcount\nR1\t42.5\t1\n")
    with pytest.raises(ValueError, match="integer"):
        load_tes_tsv(tes)
    tes.write_text("transcript_id\ttes\tcount\nR1\t42\tinf\n")
    with pytest.raises(ValueError, match="finite"):
        load_tes_tsv(tes)
    tes.write_text("transcript_id\ttes\tcount\nR1\t42\t0\nR1\t43\t0\n")
    with pytest.raises(ValueError, match="positive total"):
        load_tes_tsv(tes)


def test_conflicting_transcript_id_is_rejected(tmp_path):
    """One transcript identifier cannot alias conflicting structures or genes."""
    catalog = tmp_path / "catalog.json"
    catalog.write_text(json.dumps({"genes": {
        "G1": {"gene_id":"G1", "gene_name":"A", "chrom":"chr1", "strand":"+", "ref":[{"id":"same", "exons":[[0,10]]}], "lr":[]},
        "G2": {"gene_id":"G2", "gene_name":"B", "chrom":"chr1", "strand":"+", "ref":[{"id":"same", "exons":[[20,30]]}], "lr":[]}}}))
    with pytest.raises(ValueError, match="same.*conflicting"):
        load_regionquant_catalog(catalog)


def test_pigeon_filters_are_counted_and_exact_reference_gene_id_is_retained(tmp_path):
    """Raw classifications apply conservative filters while filtered-lite omissions remain usable."""
    gtf = tmp_path / "ref.gtf"
    gtf.write_text('chr1\tr\texon\t1\t10\t.\t+\t.\tgene_id "ENSG1.4"; gene_name "DUP"; transcript_id "R";\n'
                   'chr2\tr\texon\t1\t10\t.\t+\t.\tgene_id "ENSG2.7"; gene_name "DUP"; transcript_id "S";\n')
    gff = tmp_path / "p.gff"
    gff.write_text("".join(f'chr1\tp\texon\t{20+i*20}\t{29+i*20}\t.\t+\t.\tParent=P{i}\n' for i in range(5)))
    cls = tmp_path / "c.tsv"
    cls.write_text("isoform\tchrom\tstrand\tassociated_gene\tstructural_category\tfl_assoc\tRTS_stage\tall_canonical\n"
                   "P0\tchr1\t+\tENSG1.9\tnovel_not_in_catalog\t2\tFALSE\tcanonical\n"
                   "P1\tchr1\t+\tENSG1.9\tnovel_not_in_catalog\t2\tTRUE\tcanonical\n"
                   "P2\tchr1\t+\tENSG1.9\tnovel_not_in_catalog\t2\tFALSE\tnon_canonical\n"
                   "P3\tchr1\t+\tENSG1.9\tfusion\t2\tFALSE\tcanonical\n"
                   "P4\tchr1\t+\tNOVEL\tintergenic\t1\tFALSE\tcanonical\n")
    tx, qc = load_transcript_catalog(gtf, gff, cls)
    assert next(x for x in tx if x.id == "R").gene_id == "ENSG1.4"
    assert next(x for x in tx if x.id == "P0").gene_id == "ENSG1.4"
    assert next(x for x in tx if x.id == "P4").gene_id == "NOVEL"
    assert qc["excluded"] == {"rts": 1, "noncanonical": 1, "unsupported_category": 1}
    assert "DUP" not in qc["gene_aliases"]
