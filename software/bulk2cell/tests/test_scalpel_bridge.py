"""Verify the external SCALPEL bridge preserves identities and GTF coordinates."""
import importlib.util
from pathlib import Path
import pytest
import pysam


def load_bridge(name="scalpel_bridge"):
    """Load the standalone benchmark runner as an importable test module."""
    path = Path(__file__).parents[1]/"benchmark/scalpel/run.py"
    spec = importlib.util.spec_from_file_location(name, path)
    bridge = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bridge)
    return bridge


def test_gtf_exon_identity_and_closed_coordinates(tmp_path):
    """Standard GTF exons pass through without an off-by-one coordinate change."""
    path = Path(__file__).parents[1]/"benchmark/scalpel/run.py"
    spec = importlib.util.spec_from_file_location("scalpel_bridge", path)
    bridge = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bridge)
    gtf = tmp_path/"input.gtf"
    gtf.write_text('chr1\ttest\texon\t101\t200\t.\t-\t.\tgene_id "ENSG_1"; transcript_id "PB_1.1";\n')
    assert list(bridge.exon_rows(gtf)) == [("chr1",101,200,"-","ENSG_1","PB_1.1")]


def test_validate_prepared_reads_checks_inputs_and_bed_hashes(tmp_path):
    """A completed shared bundle is accepted only when every binding matches."""
    import json
    bridge = load_bridge("scalpel_prepared_reads")
    bundle = tmp_path/"shared"; (bundle/"bed").mkdir(parents=True)
    bam = tmp_path/"reads.bam"; bam.write_bytes(b"bam")
    barcodes = tmp_path/"barcodes.tsv"; barcodes.write_text("cell\n")
    bed = bundle/"bed/chr1.bed"; bed.write_text("chr1\t1\t2\t+\tr/1\tCB:Z:cell::UB:Z:u\n")
    def record(path):
        item = bridge.file_fingerprint(path)
        return {"path": item["path"], "bytes": item["size"], "sha256": item["sha256"]}
    certificate = {"schema":"bulk2cell.shared_reads.v1", "completed":True,
        "cross_cb_read_ids":0, "projection_records":1,
        "selection_rules":["mapped", "primary", "not supplementary", "not QC-fail",
                           "CB in exact 10,250-cell whitelist", "nonempty UB"],
        "read_id_rule":"query_name before first /", "bam":record(bam),
        "canonical_sorted_barcodes":record(barcodes),
        "scalpel_consumer_source":record(Path(__file__).parents[1]/"benchmark/scalpel/run.py"),
        "bed_files":{"chr1.bed":record(bed)}}
    (bundle/"certificate.json").write_text(json.dumps(certificate))
    result = bridge.validate_prepared_reads(bundle, bam, barcodes, {"chr1"})
    assert result["bed_paths"]["chr1"] == bed
    bed.write_text("changed")
    with pytest.raises(ValueError, match="fingerprint"):
        bridge.validate_prepared_reads(bundle, bam, barcodes, {"chr1"})


def test_prepare_chromosome_shard_preserves_complete_cross_gene_umi_family(tmp_path):
    """All records in a chromosome-level CB::UB family remain visible upstream."""
    bridge = load_bridge("scalpel_chromosome")
    gtf = tmp_path/"input.gtf"
    gtf.write_text(
        'chr1\tt\texon\t101\t130\t.\t+\t.\tgene_id "g1"; transcript_id "t1";\n'
        'chr1\tt\texon\t151\t180\t.\t+\t.\tgene_id "g2"; transcript_id "t2";\n')
    quant = tmp_path/"quant.sf"
    quant.write_text("Name\tLength\tEffectiveLength\tTPM\tNumReads\nt1\t30\t20\t1\t1\nt2\t30\t20\t1\t1\n")
    barcodes = tmp_path/"barcodes.tsv"; barcodes.write_text("cell\n")
    bam = tmp_path/"reads.bam"
    header = {"HD":{"VN":"1.6"}, "SQ":[{"SN":"chr1","LN":1000}]}
    with pysam.AlignmentFile(bam, "wb", header=header) as out_bam:
        for name, start, cigar, length in (
                ("exonic", 100, [(0,30),(3,20),(0,30)], 60),
                ("intronic_collision", 400, [(0,20)], 20)):
            read = pysam.AlignedSegment(); read.query_name=name
            read.query_sequence="A"*length; read.flag=0; read.reference_id=0
            read.reference_start=start; read.cigar=cigar
            read.query_qualities=pysam.qualitystring_to_array("I"*length)
            read.set_tag("CB","cell"); read.set_tag("UB","same_umi"); read.set_tag("GX","g1")
            out_bam.write(read)
    out = tmp_path/"out"; out.mkdir()
    shards, qc = bridge.prepare(gtf, quant, bam, barcodes, out)
    assert len(shards) == 1
    bed = (out/shards[0]/"reads.bed").read_text().splitlines()
    assert len(bed) == 3
    assert {line.split("\t")[5] for line in bed} == {"CB:Z:cell::UB:Z:same_umi"}
    assert qc["selected_alignments"] == 2


def test_prepare_rejects_gene_split_across_chromosomes(tmp_path):
    """Per-chromosome EM is valid only when each gene has one chromosome."""
    bridge = load_bridge("scalpel_gene_chromosome")
    gtf = tmp_path/"input.gtf"
    gtf.write_text(
        'chr1\tt\texon\t1\t10\t.\t+\t.\tgene_id "g"; transcript_id "t1";\n'
        'chr2\tt\texon\t1\t10\t.\t+\t.\tgene_id "g"; transcript_id "t2";\n')
    with pytest.raises(ValueError, match="multiple chromosomes"):
        bridge.prepare(gtf, tmp_path/"missing.quant", tmp_path/"missing.bam",
                       tmp_path/"missing.barcodes", tmp_path)


def test_assign_cb_batches_keeps_cells_whole_and_balances_records():
    """Whole cells are assigned once while record-heavy cells are balanced."""
    bridge = load_bridge("scalpel_cb_plan")
    plan = bridge.assign_cb_batches({"A": 10, "B": 7, "C": 2, "D": 1}, 2)
    assert set(plan) == {"A", "B", "C", "D"}
    assert plan["A"] != plan["B"]
    assert len(set(plan.values())) == 2


def test_repartition_by_cb_preserves_whole_families_and_full_annotation(tmp_path):
    """Cell batches retain complete families and link the full chromosome GTF."""
    bridge = load_bridge("scalpel_cb_repartition")
    chrom = tmp_path/"chrom_000"; chrom.mkdir()
    (chrom/"annotation.tsv").write_text("full chromosome annotation\n")
    (chrom/"reads.bed").write_text(
        "chr1\t1\t2\t+\tr1/1\tCB:Z:A::UB:Z:u\n"
        "chr1\t3\t4\t+\tr1/2\tCB:Z:A::UB:Z:u\n"
        "chr1\t5\t6\t+\tr2/1\tCB:Z:B::UB:Z:v\n")
    shards, qc = bridge.repartition_by_cb(
        tmp_path, ["chrom_000"], {"chrom_000": "chr1"}, 1,
        {"unique_read_id_cb_pairs": 2, "cross_cb_read_ids": 0})
    assert len(shards) == 2
    contents = [(tmp_path/shard/"reads.bed").read_text() for shard in shards]
    assert sum("CB:Z:A::UB:Z:u" in content for content in contents) == 1
    assert sum(content.count("CB:Z:A::UB:Z:u") for content in contents) == 2
    assert all((tmp_path/shard/"annotation.tsv").resolve() ==
               (chrom/"annotation.tsv").resolve() for shard in shards)
    assert qc["cross_cb_readids"] == 0
    assert qc["prepared_cross_cb_certificate"] == 1
    assert not (tmp_path/"readid_cb.tsv").exists()
    assert not (chrom/"reads.bed").exists()


def test_upstream_read_id_uses_text_before_first_slash():
    """The gate matches tidyr separate semantics for embedded slash names."""
    bridge = load_bridge("scalpel_readid_semantics")
    assert bridge.upstream_read_id("instrument/read/3") == "instrument"


def test_readid_cb_gate_rejects_cross_cell_identity(tmp_path):
    """A read ID reused by different cells fails before upstream filtering."""
    bridge = load_bridge("scalpel_readid_gate")
    pairs = tmp_path/"pairs.tsv"
    pairs.write_text("r2\tB\nr1\tA\nr1\tB\nr2\tB\n")
    with pytest.raises(ValueError, match="multiple cell barcodes"):
        bridge.validate_readid_cb_pairs(pairs, tmp_path)


def test_mapping_bridges_only_instrument_and_preserve_global_noop(tmp_path):
    """Generated bridges record candidates and conditionally preserve the native bug."""
    bridge = load_bridge("scalpel_mapping_bridge")
    source = tmp_path/"mapping.R"
    source.write_text("before\n#---- concordance unspliced/spliced fragments\ntrs.todel = thing\nunspliceds = dplyr::filter(unspliceds, !(ftrs %in% trs.todel))\nafter\n")
    instrumented, noop = bridge.write_mapping_bridges(source, tmp_path)
    assert "fwrite(trs.todel" in instrumented.read_text()
    assert "dplyr::filter(unspliceds" in instrumented.read_text()
    assert "unspliceds = unspliceds" in noop.read_text()
    assert "dplyr::filter(unspliceds" not in noop.read_text()


def test_concordance_plan_reruns_only_singletons_when_global_count_exceeds_one(tmp_path):
    """The bridge recreates `%in% data.table` cardinality semantics exactly."""
    bridge = load_bridge("scalpel_concordance_plan")
    paths = []
    for name, text in (("a", "f1\n"), ("b", "f2\n"), ("c", "")):
        path = tmp_path/f"{name}.tsv"; path.write_text(text); paths.append(path)
    result = bridge.concordance_correction_plan(paths, tmp_path/"global.tsv", tmp_path)
    assert result["global_candidates"] == 2
    assert result["rerun_indexes"] == [0, 1]


def test_reconcile_observed_pairs_uses_global_transcript_count(tmp_path):
    """Genes unique within batches are rejected when globally multi-model."""
    bridge = load_bridge("scalpel_global_pairs")
    first = tmp_path/"a.tsv"; second = tmp_path/"b.tsv"
    first.write_text("g1\tt1\ng2\tt2\n")
    second.write_text("g1\tt3\ng2\tt2\n")
    genes = tmp_path/"unique_genes.tsv"
    result = bridge.reconcile_observed_pairs([first, second], genes, tmp_path)
    assert result == {"observed_pairs": 3, "unique_genes": 1}
    assert genes.read_text() == "g2\n"


def test_batch_r_helpers_reconcile_nonempty_original_unique_rows(tmp_path):
    """R helpers extract observed pairs and apply a global gene allowlist."""
    import subprocess
    root = Path(__file__).parents[1]/"benchmark/scalpel"
    rscript = Path(__file__).parents[1]/".tools/scalpel/env/bin/Rscript"
    if not rscript.exists():
        pytest.skip("pinned SCALPEL R environment unavailable")
    filtered = tmp_path/"filtered.rds"
    subprocess.run([rscript, "-e", (
        'saveRDS(data.frame(seqnames.rd="chr1",start.rd=1,end.rd=2,strand.rd="+",'
        'dist_END=3,frag.id=c("c::u1","c::u2"),start=1,end=9,gene_name=c("g1","g2"),'
        'transcript_name=c("t1","t2"),bulk_TPMperc=.5),commandArgs(TRUE)[1])'), filtered], check=True)
    pairs = tmp_path/"pairs.tsv"
    subprocess.run([rscript, root/"batch_pairs.R", filtered, pairs], check=True)
    assert set(pairs.read_text().splitlines()) == {"g1\tt1", "g2\tt2"}
    genes = tmp_path/"genes.tsv"; genes.write_text("g2\n")
    unique = tmp_path/"unique.reads"
    subprocess.run([rscript, root/"batch_unique.R", filtered, genes, unique], check=True)
    assert "g2\tt2" in unique.read_text()
    assert "g1\tt1" not in unique.read_text()


def test_empty_filtered_batch_is_explicitly_ineligible_for_fragment_r(tmp_path):
    """A successful zero-row filtered batch bypasses fread and EM explicitly."""
    bridge = load_bridge("scalpel_empty_filtered")
    pairs = tmp_path/"observed_pairs.tsv"; pairs.write_text("")
    assert bridge.filtered_batch_has_evidence(pairs) is False
    pairs.write_text("g\tt\n")
    assert bridge.filtered_batch_has_evidence(pairs) is True


def test_diagnostic_sink_is_explicit_devnull_symlink(tmp_path):
    """The unused upstream distance diagnostic is discarded with provenance."""
    bridge = load_bridge("scalpel_diagnostic_sink")
    bridge.install_diagnostic_sink(tmp_path)
    diagnostic = tmp_path/"read_distance_distribution_on_transcriptomic_scope.txt"
    assert diagnostic.is_symlink()
    assert diagnostic.resolve() == Path("/dev/null")
    assert "discard" in (tmp_path/"diagnostic_provenance.json").read_text()


def test_resume_rejects_changed_parameters_and_missing_completed_outputs(tmp_path):
    """Resume state is bound to parameters and proves every declared output exists."""
    bridge = load_bridge("scalpel_resume")
    source = tmp_path/"source"; source.write_text("input")
    output = tmp_path/"artifact"; output.write_text("done")
    manifest = tmp_path/"stage.json"
    identity = bridge.stage_identity("stage", {"source": source}, {"binsize": 20})
    bridge.write_stage_manifest(manifest, identity, [output], {"seconds": 1.0, "peak_rss_kib": 10})
    assert bridge.completed_stage(manifest, identity) is not None
    changed = bridge.stage_identity("stage", {"source": source}, {"binsize": 21})
    with pytest.raises(ValueError, match="does not match"):
        bridge.completed_stage(manifest, changed)
    output.unlink()
    with pytest.raises(FileNotFoundError, match="declared output"):
        bridge.completed_stage(manifest, identity)


def test_run_lock_rejects_concurrent_writer(tmp_path):
    """Only one process may mutate a SCALPEL run directory at a time."""
    bridge = load_bridge("scalpel_lock")
    first = bridge.acquire_run_lock(tmp_path)
    try:
        with pytest.raises(RuntimeError, match="already active"):
            bridge.acquire_run_lock(tmp_path)
    finally:
        first.close()


def test_sparse_export_preserves_zero_cells_and_collapsed_members(tmp_path):
    """Collapsed memberships and zero-cell alignment survive matrix conversion."""
    import json
    from scipy.io import mmread
    path = Path(__file__).parents[1]/"benchmark/scalpel/export.py"
    spec = importlib.util.spec_from_file_location("scalpel_export", path)
    bridge = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bridge)
    (tmp_path/"chr1").mkdir()
    (tmp_path/"aliases.json").write_text(json.dumps({"G1":"gene_a", "T1":"tx_a", "T2":"tx_b"}))
    (tmp_path/"summary.json").write_text("{}")
    (tmp_path/"chr1/exons.tsv").write_text("transcript_name\tgene_name\tcollapsed\nT1\tG1\tT1_T2\n")
    (tmp_path/"estimates.tsv").write_text("bc\tgene_name\ttranscript_name\testimated_count\nCELL1\tgene_a\ttx_a\t2.5\n")
    barcodes = tmp_path/"selected.tsv"
    barcodes.write_text("CELL2\nCELL1\n")
    result = bridge.export(tmp_path, barcodes)
    assert result["cells"] == 2
    assert mmread(tmp_path/"isoform_counts.mtx.gz").toarray().tolist() == [[2.5], [0.0]]
    assert "tx_a\tgene_a\ttx_b" in (tmp_path/"membership.tsv").read_text()


def test_sparse_export_refuses_incomplete_run_without_writing(tmp_path):
    """An interrupted upstream run cannot publish benchmark matrices."""
    path = Path(__file__).parents[1]/"benchmark/scalpel/export.py"
    spec = importlib.util.spec_from_file_location("scalpel_export_incomplete", path)
    bridge = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bridge)
    barcodes = tmp_path/"selected.tsv"
    barcodes.write_text("CELL1\n")

    with pytest.raises(FileNotFoundError, match="summary.json"):
        bridge.export(tmp_path, barcodes)

    assert not any((tmp_path/name).exists() for name in (
        "isoform_counts.mtx.gz", "gene_counts.mtx.gz", "barcodes.tsv",
        "isoforms.tsv", "genes.tsv", "membership.tsv", "matrix_summary.json",
    ))


def test_resume_rejects_changed_completed_output_bytes(tmp_path):
    """A present but altered output cannot satisfy a resumable stage manifest."""
    bridge = load_bridge("scalpel_resume_output_hash")
    source = tmp_path/"source"; source.write_text("input")
    output = tmp_path/"artifact"; output.write_text("exact")
    identity = bridge.stage_identity("stage", {"source": source}, {"binsize": 20})
    manifest = tmp_path/"stage.complete.json"
    bridge.write_stage_manifest(manifest, identity, [output],
                                {"seconds": 1.0, "peak_rss_kib": 10})
    output.write_text("tampered")
    with pytest.raises(ValueError, match="fingerprint changed"):
        bridge.completed_stage(manifest, identity)


def test_bounded_workers_fail_fast_and_cancel_queued_commands(tmp_path):
    """The first failed command terminates its sibling and prevents queued launches."""
    import time
    bridge = load_bridge("scalpel_worker_failure")
    started_markers = []

    def task(index):
        """Run one owned subprocess with one deliberate early failure."""
        if index == 0:
            command = ["/bin/sh", "-c", "sleep 5"]
        elif index == 1:
            command = ["/bin/sh", "-c", "sleep 0.2; exit 7"]
        else:
            marker = tmp_path/f"{index}.started"
            started_markers.append(marker)
            command = ["/bin/sh", "-c", f"touch {marker}; sleep 5"]
        return bridge.run_command(command, tmp_path, tmp_path/f"{index}.log")

    started = time.monotonic()
    with pytest.raises(Exception):
        bridge.bounded_map(task, range(8), 2)
    assert time.monotonic() - started < 2
    assert not any(path.exists() for path in started_markers)
    assert not bridge._ACTIVE_CHILDREN


def test_invocation_ledger_accumulates_completed_and_flags_gap(tmp_path):
    """Completed attempts accumulate while an unclosed attempt remains explicit."""
    import json
    bridge = load_bridge("scalpel_invocation_ledger")
    first, prior, gaps = bridge.start_invocation_ledger(tmp_path, {"resume": False})
    assert (prior, gaps) == (0, 0)
    bridge.finish_invocation_ledger(tmp_path, first, 12.5)
    bridge.start_invocation_ledger(tmp_path, {"resume": True})
    third, prior, gaps = bridge.start_invocation_ledger(tmp_path, {"resume": True})
    assert (prior, gaps) == (12.5, 1)
    bridge.finish_invocation_ledger(tmp_path, third, 3.5)
    rows = json.loads((tmp_path/"invocations.json").read_text())["invocations"]
    assert [row["status"] for row in rows] == ["complete", "running", "complete"]
    assert sum(row.get("wall_seconds", 0) for row in rows) == 16


def test_adoption_lineage_normalizes_bytes_to_fingerprint_size():
    """Certificate lineage bytes match the equivalent stage-fingerprint size."""
    path = Path(__file__).parents[1]/"benchmark/scalpel/adopt_prepared.py"
    spec = importlib.util.spec_from_file_location("scalpel_adoption", path)
    adoption = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adoption)
    lineage = {"path": "/x", "bytes": 4, "sha256": "abc"}
    fingerprint = {"path": "/x", "size": 4, "sha256": "abc"}
    assert adoption.lineage_matches_fingerprint(lineage, fingerprint)
    assert not adoption.lineage_matches_fingerprint(
        {**lineage, "bytes": 5}, fingerprint)


def test_adoption_lineage_follows_exact_chain_and_rejects_cycle(tmp_path):
    """Every ancestry hop verifies its certificate and cycles are rejected."""
    import json
    path = Path(__file__).parents[1]/"benchmark/scalpel/adopt_prepared.py"
    spec = importlib.util.spec_from_file_location("scalpel_adoption_chain", path)
    adoption = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(adoption)
    original = tmp_path/"original"; parent = tmp_path/"parent"; child = tmp_path/"child"
    for directory, content in ((original, "old"), (parent, "parent"), (child, "child")):
        directory.mkdir()
        (directory/"certificate.json").write_text(content)

    def lineage_record(certificate):
        """Translate a local file fingerprint to certificate lineage schema."""
        record = adoption.sha256_record(certificate)
        return {"path": record["path"], "bytes": record["size"],
                "sha256": record["sha256"]}

    old = adoption.sha256_record(original/"certificate.json")
    (parent/"lineage.json").write_text(json.dumps({
        "source_bundle": str(original),
        "source_certificate": lineage_record(original/"certificate.json")}))
    (child/"lineage.json").write_text(json.dumps({
        "source_bundle": str(parent),
        "source_certificate": lineage_record(parent/"certificate.json")}))
    evidence = adoption.verified_lineage_chain(child/"lineage.json", old)
    assert evidence and evidence[-1] == old

    unrelated = tmp_path/"unrelated"; unrelated.mkdir()
    (unrelated/"certificate.json").write_text("unrelated")
    (unrelated/"lineage.json").write_text((parent/"lineage.json").read_text())
    forged = json.loads((child/"lineage.json").read_text())
    forged["source_bundle"] = str(unrelated)
    (child/"lineage.json").write_text(json.dumps(forged))
    with pytest.raises(ValueError, match="hop is invalid"):
        adoption.verified_lineage_chain(child/"lineage.json", old)

    (child/"lineage.json").write_text(json.dumps({
        "source_bundle": str(parent),
        "source_certificate": lineage_record(parent/"certificate.json")}))
    (parent/"lineage.json").write_text(json.dumps({
        "source_bundle": str(child),
        "source_certificate": lineage_record(child/"certificate.json")}))
    with pytest.raises(ValueError, match="cycle"):
        adoption.verified_lineage_chain(
            child/"lineage.json", {"path": "/missing", "size": 0, "sha256": "x"})
