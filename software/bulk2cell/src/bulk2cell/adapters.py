"""Adapters for annotation, abundance, TES, barcode, and indexed 10x BAM inputs.

All genomic intervals are zero-based half-open. BAM extraction requires a called
CB, corrected UB, exactly one matching GX, and a primary sense alignment. Missing
GX records are counted and excluded; no novel gene assignment is inferred.
Raw Pigeon classifications retain FSM, ISM, NIC, NNC, intergenic, genic, and
genic-intron categories, while RTS, explicitly noncanonical, fusion, antisense,
and other categories are excluded. Missing RTS/canonical columns are accepted so
an already filtered-lite classification remains readable.
"""
from __future__ import annotations
import csv, gzip, json, math, re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable, Mapping
from .models import Molecule, Transcript

_ATTR = re.compile(r'([^\s;=]+)\s+"([^"]*)"')
_ENSG = re.compile(r"^(ENS[A-Z]*G\d+)\.\d+(_PAR_Y)?$")

def _open(path):
    """Open plain or gzip-compressed text."""
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "rt")

def normalize_gene_id(identifier: str) -> str:
    """Remove numeric Ensembl versions without damaging dotted gene symbols."""
    value = identifier.strip(); match = _ENSG.fullmatch(value)
    return match.group(1) + (match.group(2) or "") if match else value

def _attrs(text):
    """Parse quoted GTF and simple key-value GFF3 attributes."""
    result = dict(_ATTR.findall(text))
    for item in text.strip().rstrip(";").split(";"):
        if "=" in item:
            key, value = item.strip().split("=", 1); result.setdefault(key, value)
    return result

def _exons(path):
    """Yield exon rows converted from one-based closed coordinates."""
    with _open(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"): continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9: raise ValueError(f"{path}:{number}: expected nine columns")
            if fields[2].lower() != "exon": continue
            attrs = _attrs(fields[8]); tid = attrs.get("transcript_id") or attrs.get("Parent")
            if tid and tid.startswith("transcript:"): tid = tid[11:]
            if not tid: raise ValueError(f"{path}:{number}: exon lacks transcript identifier")
            start, end = int(fields[3])-1, int(fields[4])
            if start < 0 or end <= start or fields[6] not in "+-": raise ValueError(f"{path}:{number}: invalid exon")
            yield tid, fields[0], fields[6], start, end, attrs

def _dedup(records):
    """Deduplicate exact structures, rejecting reused IDs with conflicting identity."""
    groups = defaultdict(list)
    identities = {}
    for row in records:
        identity = (row["gene_id"], row["chrom"], row["strand"], row["exons"])
        previous = identities.setdefault(row["id"], identity)
        if previous != identity:
            raise ValueError(f"transcript ID {row['id']} has conflicting gene, locus, or structure")
        groups[identity].append(row)
    output, aliases, removed = [], {}, 0
    for (gid, chrom, strand, exons), rows in groups.items():
        primary = next((r for r in rows if r["source"] == "reference"), rows[0]); cid = primary["id"]
        for row in rows: aliases[row["id"]] = cid
        count = sum(r["lr_count"] for r in rows if r["source"] != "reference")
        source = "+".join(sorted({r["source"] for r in rows}))
        output.append(Transcript(cid, gid, chrom, strand, exons, count, source)); removed += len(rows)-1
    output.sort(key=lambda x: (x.chrom, x.exons[0][0], x.gene_id, x.id))
    return output, aliases, removed

def load_transcript_catalog(reference_gtf, pigeon_gff=None, pigeon_classification=None):
    """Load a filtered reference/Pigeon union, retaining novel genes and aliases.

    Accepted Pigeon categories are FSM, ISM, NIC, NNC, intergenic, genic, and
    genic-intron. Rows explicitly marked RTS or noncanonical are excluded. Omitted
    filter columns are permitted for Pigeon filtered-lite inputs.
    """
    if (pigeon_gff is None) != (pigeon_classification is None): raise ValueError("both Pigeon inputs are required")
    records, models, gene_aliases, excluded = [], {}, defaultdict(set), Counter()
    for tid, chrom, strand, start, end, attrs in _exons(reference_gtf):
        if "gene_id" not in attrs: raise ValueError(f"reference transcript {tid} lacks gene_id")
        gid = attrs["gene_id"].strip(); key = (tid, gid)
        row = models.setdefault(key, dict(id=tid,gene_id=gid,chrom=chrom,strand=strand,exons=[],lr_count=0.,source="reference"))
        if (row["chrom"],row["strand"]) != (chrom,strand): raise ValueError(f"transcript {tid} spans loci")
        row["exons"].append((start,end))
        for alias in (gid, normalize_gene_id(gid), attrs.get("gene_name",gid)): gene_aliases[(alias,chrom,strand)].add(gid)
    if pigeon_gff is not None:
        structures = {}
        for tid, chrom, strand, start, end, attrs in _exons(pigeon_gff):
            row = structures.setdefault(tid, dict(chrom=chrom,strand=strand,exons=[]))
            if (row["chrom"],row["strand"]) != (chrom,strand): row["conflict"] = True
            row["exons"].append((start,end))
        with _open(pigeon_classification) as handle:
            seen_classifications = set()
            for row in csv.DictReader(handle, delimiter="\t"):
                tid=row.get("isoform",""); structure=structures.get(tid)
                if tid in seen_classifications: raise ValueError(f"duplicate Pigeon classification row: {tid}")
                seen_classifications.add(tid)
                category = row.get("structural_category", "")
                accepted = {"full-splice_match", "incomplete-splice_match", "novel_in_catalog", "novel_not_in_catalog", "intergenic", "genic", "genic_intron"}
                if category not in accepted: excluded["unsupported_category"] += 1; continue
                if row.get("RTS_stage", "").upper() == "TRUE": excluded["rts"] += 1; continue
                if row.get("all_canonical", "") not in ("", "canonical", "NA"): excluded["noncanonical"] += 1; continue
                if not structure or structure.get("conflict"): excluded["missing_or_conflicting_structure"]+=1; continue
                chrom=row.get("chrom") or structure["chrom"]; strand=row.get("strand") or structure["strand"]
                if (chrom,strand)!=(structure["chrom"],structure["strand"]): excluded["classification_locus_conflict"]+=1; continue
                associated=row.get("associated_gene","").strip(); normalized=normalize_gene_id(associated)
                hits=gene_aliases[(associated,chrom,strand)] | gene_aliases[(normalized,chrom,strand)]
                if len(hits)>1: excluded["ambiguous_gene_mapping"]+=1; continue
                gid=next(iter(hits),normalized)
                if not gid: excluded["missing_gene_mapping"]+=1; continue
                try: count=float(row.get("fl_assoc",0) or 0)
                except ValueError as exc: raise ValueError(f"invalid fl_assoc for {tid}") from exc
                if count<0: raise ValueError(f"negative fl_assoc for {tid}")
                models[(tid,gid,"pigeon")]=dict(id=tid,gene_id=gid,chrom=chrom,strand=strand,exons=structure["exons"],lr_count=count,source="pigeon")
    for row in models.values():
        row["exons"]=tuple(sorted(set(row["exons"])))
        if row["exons"]: records.append(row)
    tx, aliases, removed = _dedup(records)
    global_aliases = defaultdict(set)
    for (alias, _, _), gids in gene_aliases.items(): global_aliases[alias].update(gids)
    gene_names = {alias: next(iter(gids)) for alias, gids in global_aliases.items() if len(gids)==1}
    return tx, dict(aliases=aliases,gene_aliases=gene_names,excluded=dict(excluded),deduplicated_transcripts=removed,transcripts=len(tx))

def load_regionquant_catalog(path):
    """Adapt a regionquant JSON catalog without importing regionquant."""
    with _open(path) as handle: payload=json.load(handle)
    if not isinstance(payload.get("genes"),dict): raise ValueError("catalog requires genes")
    rows=[]; gene_names={}
    for key,gene in payload["genes"].items():
        gid=str(gene.get("gene_id",key)); gene_names.setdefault(gene.get("gene_name",gid), set()).add(gid); gene_names.setdefault(key, set()).add(gid)
        for field,source in (("ref","reference"),("lr","pigeon")):
            for model in gene.get(field,[]):
                rows.append(dict(id=model["id"],gene_id=gid,chrom=gene["chrom"],strand=gene["strand"],exons=tuple(tuple(map(int,x)) for x in model["exons"]),lr_count=float(model.get("fl",0)) if field=="lr" else 0.,source=source))
    tx,aliases,removed=_dedup(rows)
    unique_names = {alias: next(iter(gids)) for alias, gids in gene_names.items() if len(gids) == 1}
    return tx,dict(aliases=aliases,gene_aliases=unique_names,deduplicated_transcripts=removed,transcripts=len(tx))

def load_salmon_quant(path, transcripts: Iterable[Transcript], aliases: Mapping[str,str]|None=None):
    """Aggregate Salmon TPM by aliases, falling back explicitly to NumReads."""
    ids={x.id for x in transcripts}; amap=dict(aliases or {x:x for x in ids}); values={x:0. for x in ids}; qc=Counter()
    with _open(path) as handle:
        reader=csv.DictReader(handle,delimiter="\t")
        fields=set(reader.fieldnames or ())
        if "Name" not in fields or not ({"TPM", "NumReads"} & fields):
            raise ValueError("quant.sf requires Name and TPM or NumReads")
        measure = "TPM" if "TPM" in fields else "NumReads"
        seen=set()
        for row in reader:
            name=row["Name"]
            if name in seen: raise ValueError(f"duplicate transcript row: {name}")
            seen.add(name); cid=amap.get(name)
            try: value=float(row[measure])
            except (TypeError, ValueError) as exc: raise ValueError(f"invalid {measure} for {name}") from exc
            if not math.isfinite(value): raise ValueError(f"{measure} must be finite for {name}")
            if value<0: raise ValueError(f"negative {measure} for {name}")
            if cid not in values: qc["unmatched_rows"]+=1; continue
            values[cid]+=value; qc["matched_rows"]+=1
    qc["measure"] = "TPM" if measure == "TPM" else "NumReads_fallback"
    return values,dict(qc)

def load_tes_tsv(path):
    """Read integer zero-based TES boundaries and finite optional count weights."""
    values=defaultdict(list); n=0
    with _open(path) as handle:
        reader=csv.DictReader(handle,delimiter="\t")
        if not reader.fieldnames or not {"transcript_id","tes"}<=set(reader.fieldnames): raise ValueError("TES TSV requires transcript_id and tes")
        for row in reader:
            try:
                numeric=float(row["tes"]); weight=float(row["count"]) if "count" in row and row["count"] != "" else 1.0
            except (TypeError, ValueError) as exc: raise ValueError("invalid TES position or weight") from exc
            if not math.isfinite(numeric) or not math.isfinite(weight): raise ValueError("TES positions and weights must be finite")
            if not numeric.is_integer(): raise ValueError("TES position must be an integer boundary")
            pos=int(numeric)
            if pos<0 or weight<0: raise ValueError("negative TES position or weight")
            values[row["transcript_id"]].append((pos,weight)); n+=1
    for transcript_id, observations in values.items():
        if sum(weight for _, weight in observations) <= 0:
            raise ValueError(f"TES transcript {transcript_id} must have positive total weight")
    return dict(values),dict(observations=n,transcripts=len(values))

def load_barcodes(path):
    """Load filtered barcodes while preserving GEM suffixes."""
    with _open(path) as handle: return {line.strip().split("\t",1)[0] for line in handle if line.strip()}

def _shape(read):
    """Extract aligned bases and explicit CIGAR-N junctions."""
    pos=read.reference_start; blocks=[]; junctions=[]
    for op,n in read.cigartuples or ():
        if op in (0,7,8): blocks.append((pos,pos+n)); pos+=n
        elif op==2: pos+=n
        elif op==3: junctions.append((pos,pos+n)); pos+=n
        elif op not in (1,4,5,6): raise ValueError(f"unsupported CIGAR operation {op}")
    return blocks,junctions

def _union(intervals):
    """Union overlapping or adjacent half-open intervals."""
    merged=[]
    for a,b in sorted(set(intervals)):
        if merged and a<=merged[-1][1]: merged[-1][1]=max(b,merged[-1][1])
        else: merged.append([a,b])
    return tuple(map(tuple,merged))

def extract_molecules(bam_path, transcripts: Iterable[Transcript], selected_genes: Iterable[str], barcodes: set[str],
                      *, gene_normalization=None, bam_handle=None):
    """Fetch selected gene loci from an indexed BAM and union read evidence by molecule."""
    import pysam
    transcript_list = list(transcripts)
    wanted = set(selected_genes); by_gene=defaultdict(list); normalized_genes=defaultdict(set)
    for tx in transcript_list:
        normalized_genes[normalize_gene_id(tx.gene_id)].add(tx.gene_id)
        if tx.gene_id in wanted: by_gene[tx.gene_id].append(tx)
    if wanted-set(by_gene): raise ValueError("selected genes absent: "+",".join(sorted(wanted-set(by_gene))))
    evidence={}; qc=Counter()
    from contextlib import nullcontext
    if gene_normalization is not None: normalized_genes = gene_normalization
    with (nullcontext(bam_handle) if bam_handle is not None else pysam.AlignmentFile(str(bam_path),"rb")) as bam:
        if not bam.has_index(): raise ValueError("BAM must be indexed")
        for gid,models in by_gene.items():
            loci={(x.chrom,x.strand) for x in models}
            if len(loci)!=1: raise ValueError(f"gene {gid} spans loci")
            chrom,strand=next(iter(loci)); start=min(a for x in models for a,b in x.exons); end=max(b for x in models for a,b in x.exons)
            for read in bam.fetch(chrom,start,end):
                qc["fetched_reads"]+=1
                if read.is_unmapped or read.is_secondary or read.is_supplementary or read.is_qcfail: qc["nonprimary"]+=1; continue
                cb=read.get_tag("CB") if read.has_tag("CB") else ""; ub=read.get_tag("UB") if read.has_tag("UB") else ""
                if cb not in barcodes: qc["unfiltered_barcode"]+=1; continue
                if not ub: qc["missing_ub"]+=1; continue
                if not read.has_tag("GX") or not str(read.get_tag("GX")).strip(): qc["missing_gx"]+=1; continue
                genes=[x for x in str(read.get_tag("GX")).split(";") if x]
                if len(genes)!=1: qc["multigene_gx"]+=1; continue
                gx = genes[0]
                normalized_match = (normalize_gene_id(gx) == normalize_gene_id(gid)
                                    and normalized_genes[normalize_gene_id(gid)] == {gid})
                if gx != gid and not normalized_match: qc["other_gene_gx"]+=1; continue
                if ("-" if read.is_reverse else "+")!=strand: qc["wrong_strand"]+=1; continue
                blocks,junctions=_shape(read)
                if not blocks: qc["no_aligned_bases"]+=1; continue
                family=evidence.setdefault((cb,ub,gid),dict(blocks=set(),junctions=set()))
                family["blocks"].update(blocks); family["junctions"].update(junctions); qc["eligible_reads"]+=1
    molecules=[Molecule(cb,ub,gid,_union(v["blocks"]),tuple(sorted(v["junctions"]))) for (cb,ub,gid),v in sorted(evidence.items())]
    qc["molecules"]=len(molecules)
    return molecules,dict(qc)
