"""PROTOTYPE - throwaway. Writes HD829 M20_4's prepared tables in the layout decided in
"Decide the schema and file format of the prepared tables", for the Quarto feasibility prototype
(wayfinder ticket "Prototype the settled interactions in a self-contained Quarto report").

Run from the repo root:
    .venv/bin/python tools/noise_viz/prototype_quarto/write_tables.py --out <dir> [--unfiltered-bam <bam>]

Reuses the extraction code of the layout prototype (../prototype_layout/build_prototype_data.py).
"""
import argparse
import datetime
import os
import sys

import numpy as np
import pandas as pd
import pysam
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "prototype_layout"))
import build_prototype_data as b  # noqa: E402

SCHEMA_VERSION = 1


def read_bed():
    """Read the amplicon BED (b.BED) into the amplicons table.

    Returns:
        DataFrame with one row per amplicon, in BED order: name (BED column 4), chrom, start, end (1-based closed;
        the BED's 0-based start + 1, its end unchanged), pool (BED column 5, "" when absent). Lines with fewer
        than 4 columns are skipped.
    """
    rows =[f for f in (line.split() for line in open(b.BED)) if len(f) >= 4]
    return pd.DataFrame({"name": [f[3] for f in rows], "chrom": [f[0] for f in rows],
                         "start": [int(f[1]) + 1 for f in rows], "end": [int(f[2]) for f in rows],
                         "pool": [f[4] if len(f) >= 5 else "" for f in rows]})


def merged_spans(amplicons):
    """Merge overlapping or touching amplicons into the regions to pile up, so each position is counted once.

    Args:
        amplicons: the table from read_bed() (1-based closed coordinates).

    Returns:
        List of (chrom, start, end) tuples, 0-based half-open, chromosomes in BED order and sorted by start
        within a chromosome.
    """
    out = []
    for chrom, g in amplicons.groupby("chrom", sort=False):
        for s, e in sorted(zip(g.start - 1, g.end)):
            if out and out[-1][0] == chrom and s <= out[-1][2]:
                out[-1][2] = max(out[-1][2], e)
            else:
                out.append([chrom, s, e])
    return [tuple(x) for x in out]


def count_raw(bam, chrom, start, end):
    """Pile up the reads over one region and return the raw counts per strand and position.

    Same counting rules as build_prototype_data.count() (reads accepted by b.counted(); deletions counted at
    every deleted position; insertions at their left anchor), but nothing is turned into fractions: the report
    does that, so that thresholds can change without re-reading the BAM.

    Args:
        bam: open, indexed pysam.AlignmentFile.
        chrom: chromosome name.
        start, end: the region, 0-based half-open; L = end - start.

    Returns:
        (c, qsum, qn), int64 numpy arrays:
            c: shape (2, 8, L); axis 0 is the strand (0 forward, 1 reverse); axis 1 is
               [A, C, G, T, N, deletions, insertions, inserted bases (sum of insertion lengths)];
               axis 2 is the position, offset from start;
            qsum: shape (2, L), sum of the base qualities of the aligned bases;
            qn: shape (2, L), number of aligned bases, so qsum / qn is the mean base quality.
    """
    L = end - start
    c = np.zeros((2, 8, L), dtype=np.int64)
    qsum = np.zeros((2, L), dtype=np.int64)
    qn = np.zeros((2, L), dtype=np.int64)
    for read in filter(b.counted, bam.fetch(chrom, start, end)):
        s = int(read.is_reverse)
        seq = b.LUT[np.frombuffer(read.query_sequence.encode(), dtype=np.uint8)]
        qual = np.asarray(read.query_qualities, dtype=np.int64)
        rpos, qpos, mr, mq = read.reference_start, 0, [], []
        for op, ln in read.cigartuples:
            if op in (0, 7, 8):
                mr.append(np.arange(rpos, rpos + ln))
                mq.append(np.arange(qpos, qpos + ln))
                rpos += ln
                qpos += ln
            elif op == 1:
                a = rpos - 1 - start
                if 0 <= a < L:
                    c[s, 6, a] += 1
                    c[s, 7, a] += ln
                qpos += ln
            elif op == 2:
                c[s, 5, np.arange(max(rpos, start), min(rpos + ln, end)) - start] += 1
                rpos += ln
            elif op == 3:
                rpos += ln
            elif op == 4:
                qpos += ln
        if not mr:
            continue
        rp, qp = np.concatenate(mr), np.concatenate(mq)
        keep = (rp >= start) & (rp < end)
        rp, qp = rp[keep] - start, qp[keep]
        c[s, :5] += np.bincount(seq[qp] * L + rp, minlength=5 * L).reshape(5, L)
        qsum[s] += np.bincount(rp, weights=qual[qp], minlength=L).astype(np.int64)
        qn[s] += np.bincount(rp, minlength=L)
    return c, qsum, qn


def positions(bams, amplicons, spans):
    """Build the positions table (positions.parquet): raw counts for every position, dataset and strand.

    The reference bases are fetched from UCSC (network access needed).

    Args:
        bams: {dataset: open pysam.AlignmentFile}, dataset being "filtered" or "unfiltered".
        amplicons: the table from read_bed().
        spans: the regions from merged_spans().

    Returns:
        DataFrame, one row per (position, dataset, strand), with columns
            chrom, pos (1-based), dataset, strand ("+" or "-"), ref,
            A, C, G, T, N, del, ins, ins_len (raw counts as in count_raw()),
            bq_sum, bq_n (base-quality sum and number of bases),
            covering_amplicons (list of the names of the amplicons whose span contains the position).
    """
    frames = []
    for chrom, s, e in spans:
        ref = b.ucsc_sequence(chrom, s, e)
        pos = np.arange(s + 1, e + 1)
        amps = amplicons[amplicons.chrom == chrom]
        covering = [[n for n, a, z in zip(amps.name, amps.start, amps.end) if a <= p <= z] for p in pos]
        for dataset, bam in bams.items():
            c, qsum, qn = count_raw(bam, chrom, s, e)
            for si, strand in ((0, "+"), (1, "-")):
                frames.append(pd.DataFrame({
                    "chrom": chrom, "pos": pos, "dataset": dataset, "strand": strand, "ref": list(ref),
                    "A": c[si, 0], "C": c[si, 1], "G": c[si, 2], "T": c[si, 3], "N": c[si, 4],
                    "del": c[si, 5], "ins": c[si, 6], "ins_len": c[si, 7], "bq_sum": qsum[si], "bq_n": qn[si],
                    "covering_amplicons": covering}))
        print("positions", chrom, s, e, flush=True)
    return pd.concat(frames, ignore_index=True)


def reads(fastqs):
    """Build the reads table (reads.parquet): every read of each dataset's FASTQ.

    The FASTQ, not the BAM, defines the dataset: primer clipping drops some reads from the BAM.

    Args:
        fastqs: {dataset: FASTQ path}.

    Returns:
        DataFrame with columns dataset, read_id, length (bases), mean_q (mean read quality as in
        build_prototype_data.read_quality(), rounded to 3 decimals).
    """
    frames = []
    for dataset, path in fastqs.items():
        ids, lengths, qs = [], [], []
        for r in pysam.FastxFile(path):
            ids.append(r.name)
            lengths.append(len(r.sequence))
            qs.append(float(-10 * np.log10(b.ERR_P[np.frombuffer(r.quality.encode(), dtype=np.uint8) - 33].mean())))
        frames.append(pd.DataFrame({"dataset": dataset, "read_id": ids, "length": lengths, "mean_q": np.round(qs, 3)}))
    return pd.concat(frames, ignore_index=True)


def alignments(bam_paths, amplicons):
    """Build the alignments table (alignments.parquet): one row per primary alignment of each dataset's BAM.

    No MAPQ filter is applied here; the report applies it. Secondary and supplementary alignments are left out.

    Args:
        bam_paths: {dataset: BAM path}.
        amplicons: the table from read_bed().

    Returns:
        DataFrame with columns
            dataset, read_id;
            chrom, start, end: aligned span, 1-based closed (empty for unmapped reads);
            mapq;
            amplicon: the amplicon the read is an amplicon read of (b.best_amplicon() with score >=
                b.MIN_SPAN_MATCH), empty otherwise;
            span_match: the best match score (overlap / union, 4 decimals), 0 when the read overlaps no
                amplicon, empty for unmapped reads.
    """
    by_chrom ={c: list(zip(g.name, g.start - 1, g.end)) for c, g in amplicons.groupby("chrom")}
    rows = []
    for dataset, path in bam_paths.items():
        for r in pysam.AlignmentFile(path).fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary:
                continue
            if r.is_unmapped:
                rows.append((dataset, r.query_name, None, None, None, r.mapping_quality, None, None))
                continue
            best, score = b.best_amplicon(r, by_chrom.get(r.reference_name, []))
            rows.append((dataset, r.query_name, r.reference_name, r.reference_start + 1, r.reference_end,
                         r.mapping_quality, best if score >= b.MIN_SPAN_MATCH else None, round(score, 4)))
    df = pd.DataFrame(rows, columns=["dataset", "read_id", "chrom", "start", "end", "mapq", "amplicon", "span_match"])
    return df.astype({"start": "Int64", "end": "Int64"})


def joined(x, sep):
    """Turn a VCF value into one TSV field.

    Args:
        x: None, a list or tuple (e.g. INFO/CALLERS), or a single value.
        sep: separator for list values.

    Returns:
        "" for None, the items joined by sep for a list or tuple, str(x) otherwise.
    """
    if x is None:
        return ""
    if isinstance(x, (list, tuple)):
        return sep.join(str(v) for v in x)
    return str(x)


def variants(spans, known):
    """Build the variants table (variants.tsv): the sample's calls in the amplicon regions, plus known variants.

    Args:
        spans: the regions from merged_spans().
        known: output of build_prototype_data.known_variants().

    Returns:
        DataFrame, one row per variant from build_prototype_data.variants(), with columns
            chrom, pos, ref, alt, af (6 decimals; empty when not called), dp,
            filter (";"-joined), callers (","-joined),
            artifact, artifact_median, artifact_nrsd (","-joined, one value per caller),
            known_gene, known_hgvsp, known_expected_af, known_note (empty for variants that are not known),
            category (pass, filtered, known_pass, known_filtered or known_not_called).
    """
    by_key ={(k["chrom"], int(k["pos"]), k["ref"], k["alt"]): k for k in known}
    rows = []
    for chrom, s, e in spans:
        for v in b.variants(chrom, s, e, known):
            k = by_key.get((chrom, v["pos"], v["ref"], v["alt"]), {})
            art = v.get("artifact", {})
            rows.append({
                "chrom": chrom, "pos": v["pos"], "ref": v["ref"], "alt": v["alt"],
                "af": round(v["af"], 6) if v["af"] is not None else None, "dp": v["dp"],
                "filter": joined(v["filter"], ";"), "callers": joined(v["callers"], ","),
                "artifact": joined(art.get("Artifact"), ","), "artifact_median": joined(art.get("ArtifactMedian"), ","),
                "artifact_nrsd": joined(art.get("ArtifactNrSD"), ","),
                "known_gene": k.get("gene", ""), "known_hgvsp": k.get("hgvsp", ""),
                "known_expected_af": k.get("expected_af"), "known_note": k.get("note", ""), "category": v["category"]})
    return pd.DataFrame(rows)


def homopolymers(spans):
    """Build the homopolymers table (homopolymers.tsv): runs from b.HOMOPOLYMERS lying entirely inside a region.

    Args:
        spans: the regions from merged_spans().

    Returns:
        DataFrame with columns chrom, start, end (1-based closed), base, length (bp); every run in the BED
        (3 bp and more), the report applies the minimum length.
    """
    rows = []
    for line in open(b.HOMOPOLYMERS):
        c, s, e, base, ln = line.split()[:5]
        if any(c == chrom and int(s) >= a and int(e) <= z for chrom, a, z in spans):
            rows.append((c, int(s) + 1, int(e), base, int(ln)))
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "base", "length"]).drop_duplicates()


def exons(spans):
    """Build the exons table (exons.tsv) from UCSC RefSeq Select (network access needed).

    Args:
        spans: the regions from merged_spans().

    Returns:
        DataFrame with columns chrom, gene, transcript, exon (number), strand, start, end (1-based closed,
        clipped to the region), as in build_prototype_data.exons(). The chrom column is not in the approved
        schema yet.
    """
    rows = []
    for chrom, s, e in spans:
        for x in b.exons(chrom, s, e):
            rows.append((chrom, x["gene"], x["transcript"], x["number"], x["strand"], x["start"], x["end"]))
    return pd.DataFrame(rows, columns=["chrom", "gene", "transcript", "exon", "strand", "start", "end"]).drop_duplicates()


def file_info(path):
    """Describe an input file for the provenance record in sample-metadata.yaml.

    Args:
        path: file path.

    Returns:
        {"path": absolute path, "size": bytes, "mtime": last modification time, ISO 8601 to the second}.
    """
    st =os.stat(path)
    return {"path": os.path.abspath(path), "size": st.st_size,
            "mtime": datetime.datetime.fromtimestamp(st.st_mtime).isoformat(timespec="seconds")}


def main():
    """Write the prepared tables of HD829 M20_4 to <out>/<sample>/.

    Files written: amplicons.tsv, homopolymers.tsv, exons.tsv, variants.tsv, positions.parquet, reads.parquet,
    alignments.parquet and sample-metadata.yaml (schema version, inputs with size and date, parameters, reference
    sources), then prints each file's size. The filtered BAM, FASTQ, VCF and BED paths come from
    build_prototype_data; --unfiltered-bam adds the unfiltered dataset.
    """
    ap =argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", required=True, help="prepared-tables folder; the sample's tables go to <out>/<sample>/")
    ap.add_argument("--unfiltered-bam", help="sorted, indexed, primer-clipped BAM of the unfiltered reads (optional)")
    ap.add_argument("--known-variants", default=b.KNOWN_VARIANTS)
    args = ap.parse_args()
    out = os.path.join(args.out, b.SAMPLE)
    os.makedirs(out, exist_ok=True)

    amplicons = read_bed()
    spans = merged_spans(amplicons)
    bam_paths = {"filtered": b.BAM, **({"unfiltered": args.unfiltered_bam} if args.unfiltered_bam else {})}
    bams = {k: b.open_bam(p) for k, p in bam_paths.items()}
    fastqs = {k: p for k, p in b.FASTQ.items() if k in bam_paths}

    tsv = {"sep": "\t", "index": False, "na_rep": ""}
    amplicons.to_csv(os.path.join(out, "amplicons.tsv"), **tsv)
    homopolymers(spans).to_csv(os.path.join(out, "homopolymers.tsv"), **tsv)
    exons(spans).to_csv(os.path.join(out, "exons.tsv"), **tsv)
    variants(spans, b.known_variants(args.known_variants)).to_csv(os.path.join(out, "variants.tsv"), **tsv)
    positions(bams, amplicons, spans).to_parquet(os.path.join(out, "positions.parquet"), index=False)
    reads(fastqs).to_parquet(os.path.join(out, "reads.parquet"), index=False)
    alignments(bam_paths, amplicons).to_parquet(os.path.join(out, "alignments.parquet"), index=False)

    meta = {
        "schema_version": SCHEMA_VERSION, "sample": b.SAMPLE,
        "tool": "noise_viz prototype_quarto/write_tables.py (throwaway)",
        "date": datetime.datetime.now().isoformat(timespec="seconds"),
        "inputs": {"bam": {k: file_info(p) for k, p in bam_paths.items()},
                   "fastq": {k: file_info(p) for k, p in fastqs.items()},
                   "vcf": file_info(b.VCF), "bed": file_info(b.BED), "homopolymers": file_info(b.HOMOPOLYMERS),
                   "known_variants": file_info(args.known_variants)},
        "parameters": {"min_mapq": b.MIN_MAPQ, "min_span_match": b.MIN_SPAN_MATCH},
        "reference": {"build": "GRCh38 (hg38)", "sequence": "UCSC API", "exons": "UCSC ncbiRefSeqSelect"},
    }
    with open(os.path.join(out, "sample-metadata.yaml"), "w") as f:
        yaml.safe_dump(meta, f, sort_keys=False)
    for name in sorted(os.listdir(out)):
        print(f"{name:24} {os.path.getsize(os.path.join(out, name)) / 1e6:7.2f} MB")


if __name__ == "__main__":
    main()
