"""PROTOTYPE - throwaway. Builds prototype_data.js for prototype_layout.html.

Question it answers: how should the stacked per-amplicon noise view look and behave?
(wayfinder ticket "Prototype the stacked plot layout on one HD829 amplicon").

Run from the repo root:
    .venv/bin/python tools/noise_viz/prototype_layout/build_prototype_data.py [--known-variants path.csv]
then open tools/noise_viz/prototype_layout/prototype_layout.html in a browser.

Counting follows "Decide how read noise is counted at each position": primary reads with
MAPQ >= 20, substitutions per alt base, deletions at every deleted position, insertions at
the left anchor, per strand, depth = bases + deletions. Reference bases are fetched from the
UCSC hg38 API because no GRCh38 FASTA is available locally.
"""
import argparse
import csv
import json
import os
import urllib.request

import numpy as np
import pysam

HD829 = "/home/camille/ampliconthemato/HD829"
SAMPLE = "M20_4_HD829"
BAM = f"{HD829}/M20_4/Data/{SAMPLE}_T.ont_reads.filtered.aligned.sorted.soft-clipped.bam"
VCF = f"{HD829}/M20_4/Data/{SAMPLE}_T.ont_reads.snv-indels.vcf"
FASTQ = {
    "filtered": f"{HD829}/M20_4/Data/{SAMPLE}_T.ont_reads.filtered.fastq.gz",
    "unfiltered": f"{HD829}/M20_4/Data/{SAMPLE}_T.ont_reads.unfiltered.fastq.gz",
}
BED = f"{HD829}/M20_4/Bed/amplicons.bed"
HOMOPOLYMERS = "/home/camille/ampliconthemato/amplic_ont_hemato/data/homopolymers_amplicons.bed"
KNOWN_VARIANTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "known_variants", "HD829.csv")
MIN_MAPQ, MIN_DEPTH, MIN_SPAN_MATCH = 20, 100, 0.5
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "prototype_data.js")

# Byte -> base index lookup: A=0, C=1, G=2, T=3, anything else (N, IUPAC codes, lower case) = 4.
LUT = np.full(256, 4, dtype=np.int64)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i


def ucsc_sequence(chrom, start, end):
    """Fetch a stretch of the GRCh38 (hg38) reference sequence from the UCSC REST API.

    Args:
        chrom: chromosome name with the "chr" prefix, e.g. "chr17".
        start: 0-based start of the region (included).
        end: 0-based end of the region (excluded), so the result has end - start bases.

    Returns:
        The sequence as an upper-case string (UCSC soft-masks repeats in lower case; that is dropped).

    Needs network access.
    """
    url =f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start};end={end}"
    with urllib.request.urlopen(url) as r:
        return json.load(r)["dna"].upper()


def exons(chrom, start, end):
    """List the exons overlapping a region, clipped to it, from the UCSC track "ncbiRefSeqSelect" (hg38).

    RefSeq Select has one representative transcript per gene, so each gene contributes one set of exons.

    Args:
        chrom: chromosome name with the "chr" prefix.
        start, end: the region, 0-based half-open.

    Returns:
        One dict per overlapping exon:
            gene: gene symbol, e.g. "TP53";
            transcript: RefSeq accession, e.g. "NM_000546.6";
            number: exon number in the direction of transcription (on a "-" strand gene, exon 1 is the rightmost);
            strand: "+" or "-";
            start, end: 1-based closed coordinates, clipped to the region.

    Needs network access.
    """
    url =f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=ncbiRefSeqSelect;chrom={chrom};start={start};end={end}"
    with urllib.request.urlopen(url) as r:
        txs = json.load(r)["ncbiRefSeqSelect"]
    out = []
    for tx in txs:
        starts = [int(x) for x in str(tx["exonStarts"]).strip(",").split(",")]
        ends = [int(x) for x in str(tx["exonEnds"]).strip(",").split(",")]
        for i, (s, e) in enumerate(zip(starts, ends)):
            if e > start and s < end:
                number = len(starts) - i if tx["strand"] == "-" else i + 1
                out.append({"gene": tx["name2"], "transcript": tx["name"], "number": number, "strand": tx["strand"],
                            "start": max(s, start) + 1, "end": min(e, end)})
    return out


def known_variants(path):
    """Read the known variants of a reference standard (e.g. known_variants/HD829.csv).

    Args:
        path: comma, semicolon or tab delimited file with a header row. Column names are case-insensitive.
            Required: chrom, pos, ref, alt (GRCh38, 1-based pos, VCF-style alleles; a symbolic ALT such as <INS>
            is allowed). Optional: expected_af (a fraction, e.g. 0.05), gene, hgvsc, hgvsp, note.

    Returns:
        One dict per row, with every column as a stripped string (column names lower-cased), except:
            chrom: always with the "chr" prefix;
            expected_af: float, or None when empty or ".";
            label (added): display label, "<gene> <change>", where change is hgvsp, else hgvsc, else "REF>ALT"
                (a symbolic ALT shows as e.g. "INS"); without a gene, "<chrom>:<pos>" replaces it.
        Rows with an empty chrom or pos are dropped.

    Raises:
        SystemExit: when a required column is missing.
    """
    with open(path, newline="") as f:
        dialect = csv.Sniffer().sniff(f.readline(), delimiters=",;\t")
        f.seek(0)
        rows = [{k.strip().lower(): (v or "").strip() for k, v in r.items() if k} for r in csv.DictReader(f, dialect=dialect)]
    missing = {"chrom", "pos", "ref", "alt"} - set(rows[0] if rows else {})
    if missing:
        raise SystemExit(f"{path}: missing required column(s) {', '.join(sorted(missing))}")
    for r in rows:
        r["chrom"] = r["chrom"] if r["chrom"].startswith("chr") else "chr" + r["chrom"]
        r["expected_af"] = float(r["expected_af"]) if r.get("expected_af") not in (None, "", ".") else None
        # gene + protein change, else cDNA change, else the alleles; a symbolic ALT such as <INS> becomes INS
        change = next((x for x in (r.get("hgvsp"), r.get("hgvsc")) if x and x != "."), None) or (
            r["alt"].strip("<>") if r["alt"].startswith("<") else f"{r['ref']}>{r['alt']}")
        gene = r.get("gene") if r.get("gene") not in (None, "", ".") else f"{r['chrom']}:{r['pos']}"
        r["label"] = f"{gene} {change}"
    return [r for r in rows if r["chrom"] and r["pos"]]


def count(bam, chrom, start, end, ref):
    """Pile up the reads over one region and compute the read noise at each position.

    Only reads accepted by counted() are used (primary alignments with MAPQ >= MIN_MAPQ). Each read's CIGAR is
    walked once:
        - aligned bases (M, =, X) count the read's base (A, C, G, T or N) at that reference position;
        - a deletion (D) counts once at every deleted reference position;
        - an insertion (I) counts once at its left anchor, i.e. the reference position just before it, and its
          length is added to the inserted-bases total there;
        - skipped regions (N) and soft clips (S) are stepped over; hard clips (H) and padding (P) take no space.
    Forward and reverse strand reads are counted separately.

    Definitions per position:
        depth = A + C + G + T + N + deletions (insertions sit between reference bases, so they add no depth);
        substitutions = A, C, G or T bases that differ from the reference base (N is neither match nor substitution);
        total noise = (substitutions + deletions + insertions) / depth.

    Args:
        bam: open, indexed pysam.AlignmentFile.
        chrom: chromosome name.
        start, end: the region, 0-based half-open.
        ref: reference sequence of the region (end - start bases), e.g. from ucsc_sequence().

    Returns:
        A dict; each list has one value per position of the region:
            n_reads: number of reads counted;
            pos: 1-based positions; ref: reference bases;
            depth: depth over both strands; low_depth: depth < MIN_DEPTH;
            sub, del, ins, total: noise fractions of depth, rounded to 5 decimals;
            sub_by_base: {"A": [...], "C": ..., "G": ..., "T": ...}, fraction of depth showing each non-reference base;
            ins_mean_len: mean inserted length at insertion anchors (0 where there are none);
            fwd, rev: raw counts per strand: depth, sub, del, ins;
            base_quality: mean Phred quality of the aligned bases (deletions have no quality; 0 where no bases).
        Fractions divide by max(depth, 1), so positions without reads give 0 rather than NaN.
    """
    L = end - start
    ref_idx = LUT[np.frombuffer(ref.encode(), dtype=np.uint8)]
    c = np.zeros((2, 8, L), dtype=np.int64)  # strand x [A C G T N DEL INS INSLEN] x pos
    qsum = np.zeros(L)
    qn = np.zeros(L)
    n_reads = 0
    for read in bam.fetch(chrom, start, end):
        if read.is_secondary or read.is_supplementary or read.mapping_quality < MIN_MAPQ:
            continue
        n_reads += 1
        s = int(read.is_reverse)
        seq = LUT[np.frombuffer(read.query_sequence.encode(), dtype=np.uint8)]
        qual = np.asarray(read.query_qualities)
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
                d = np.arange(max(rpos, start), min(rpos + ln, end)) - start
                c[s, 5, d] += 1
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
        qsum += np.bincount(rp, weights=qual[qp], minlength=L)
        qn += np.bincount(rp, minlength=L)

    def per(strands):
        """Reduce the counts of one strand (0 forward, 1 reverse) or both ([0, 1]).

        Returns (depth, substitutions, deletions, insertions, inserted bases, A/C/G/T counts), all raw counts
        per position; the last one is a 4 x L array.
        """
        x = c[strands].sum(axis=0) if isinstance(strands, list) else c[strands]
        depth = x[:5].sum(axis=0) + x[5]
        obs = x[:4]
        match = np.where(ref_idx < 4, obs[np.clip(ref_idx, 0, 3), np.arange(L)], 0)
        sub = obs.sum(axis=0) - match
        return depth, sub, x[5], x[6], x[7], obs

    depth, sub, dele, ins, inslen, obs = per([0, 1])
    fdepth, fsub, fdel, fins, _, _ = per(0)
    rdepth, rsub, rdel, rins, _, _ = per(1)
    d = np.maximum(depth, 1)
    r3 = lambda a: [round(float(v), 5) for v in a]  # noqa: E731  (array -> list of floats rounded to 5 decimals)
    return {
        "n_reads": n_reads,
        "pos": list(range(start + 1, end + 1)),  # 1-based
        "ref": list(ref),
        "depth": depth.tolist(),
        "low_depth": (depth < MIN_DEPTH).tolist(),
        "sub": r3(sub / d), "del": r3(dele / d), "ins": r3(ins / d),
        "total": r3((sub + dele + ins) / d),
        "sub_by_base": {b: r3(np.where(ref_idx == i, 0, obs[i]) / d) for i, b in enumerate("ACGT")},
        "ins_mean_len": r3(np.where(ins > 0, inslen / np.maximum(ins, 1), 0)),
        "fwd": {"depth": fdepth.tolist(), "sub": fsub.tolist(), "del": fdel.tolist(), "ins": fins.tolist()},
        "rev": {"depth": rdepth.tolist(), "sub": rsub.tolist(), "del": rdel.tolist(), "ins": rins.tolist()},
        "base_quality": [round(float(v), 2) for v in qsum / np.maximum(qn, 1)],
    }


def variants(chrom, start, end, known):
    """List the sample's variant calls in a region (from VCF) and mark the known variants among them.

    A call matches a known variant when chrom, pos, ref and its first ALT are equal. Known variants of the region
    that match no call are appended with category "known_not_called". Only the first ALT of a multi-allelic
    record is used.

    Args:
        chrom: chromosome name.
        start, end: the region, 0-based half-open (a variant is inside when start < pos <= end).
        known: output of known_variants().

    Returns:
        One dict per variant:
            pos (1-based), ref, alt;
            af: the sample's FORMAT/AF (first value), None when missing or not called;
            dp: the sample's FORMAT/DP;
            filter: list of FILTER values, ["PASS"] when none is set ([] when not called);
            callers: INFO/CALLERS (tuple of caller names), None when not called;
            category: "pass", "filtered" (any FILTER other than PASS), "known_pass", "known_filtered" or
                "known_not_called";
            label, expected_af: from the known variant ("" and None otherwise);
            artifact: INFO Artifact, ArtifactMedian and ArtifactNrSD when present (absent for "known_not_called").
    """
    out, seen = [], set()
    vcf = pysam.VariantFile(VCF)
    for rec in vcf.fetch(chrom, start, end) if vcf.index else vcf:
        if rec.chrom != chrom or not (start < rec.pos <= end):
            continue
        smp = rec.samples[0]
        af = smp.get("AF")
        af = af[0] if isinstance(af, tuple) else af
        filt = list(rec.filter.keys()) or ["PASS"]
        key = (rec.chrom, rec.pos, rec.ref, rec.alts[0])
        k = next((kv for kv in known if (kv["chrom"], int(kv["pos"]), kv["ref"], kv["alt"]) == key), None)
        if k:
            seen.add(key)
        cat = ("known_pass" if filt == ["PASS"] else "known_filtered") if k else ("pass" if filt == ["PASS"] else "filtered")
        out.append({
            "pos": rec.pos, "ref": rec.ref, "alt": rec.alts[0], "af": float(af) if af is not None else None,
            "dp": smp.get("DP"), "filter": filt, "callers": rec.info.get("CALLERS"),
            "category": cat, "label": k["label"] if k else "",
            "expected_af": k["expected_af"] if k else None,
            "artifact": {f: rec.info.get(f) for f in ("Artifact", "ArtifactMedian", "ArtifactNrSD") if f in rec.info},
        })
    for kv in known:
        key = (kv["chrom"], int(kv["pos"]), kv["ref"], kv["alt"])
        if kv["chrom"] == chrom and start < int(kv["pos"]) <= end and key not in seen:
            out.append({
                "pos": int(kv["pos"]), "ref": kv["ref"], "alt": kv["alt"], "af": None, "dp": None,
                "filter": [], "callers": None, "category": "known_not_called",
                "label": kv["label"], "expected_af": kv["expected_af"],
            })
    return out


def homopolymers(chrom, start, end):
    """List the homopolymer runs lying entirely inside a region, from the HOMOPOLYMERS BED.

    The BED columns used are chrom, start (0-based), end, base, length. No minimum length is applied here: the
    BED holds runs of 3 bp and more, and the viewer filters by length.

    Args:
        chrom: chromosome name.
        start, end: the region, 0-based half-open.

    Returns:
        One dict per run: start, end (1-based closed), base (e.g. "A"), length (bp).
    """
    hp = []
    for line in open(HOMOPOLYMERS):
        c, s, e, base, ln = line.split()[:5]
        if c == chrom and int(s) >= start and int(e) <= end:
            hp.append({"start": int(s) + 1, "end": int(e), "base": base, "length": int(ln)})
    return hp


ERR_P = 10 ** (-np.arange(94) / 10)  # Phred Q (index) -> error probability
Q_BINS = np.arange(0, 50.5, 0.5)  # read-quality histogram bin edges: 0, 0.5, ..., 50


def q_histogram(qs):
    """Histogram of read qualities over Q_BINS.

    Args:
        qs: list of mean read qualities (Phred), e.g. from read_quality(). Must be a list, not a numpy array.

    Returns:
        dict with counts (one per bin, 100 bins of 0.5 Q from 0 to 50; values outside are not counted),
        n (number of reads) and median (rounded to 2 decimals, None when qs is empty).
    """
    h, _ =np.histogram(qs, bins=Q_BINS)
    return {"counts": h.tolist(), "n": len(qs), "median": round(float(np.median(qs)), 2) if qs else None}


def read_quality():
    """Histogram the mean quality of every read in the filtered and the unfiltered FASTQ (see FASTQ).

    A read's mean quality is the Phred value of its mean per-base error probability,
    -10 * log10(mean(10 ** (-Q / 10))), not the mean of its Q values; it is dominated by the worst bases, so it
    is lower than the plain mean. Reads every record of both files, so it takes a while on large files.

    Returns:
        dict with bins (left bin edges) and, for "filtered" and "unfiltered", the q_histogram() of that FASTQ.
    """
    out ={"bins": Q_BINS[:-1].tolist()}
    for name, path in FASTQ.items():
        qs = [-10 * np.log10(ERR_P[np.frombuffer(r.quality.encode(), dtype=np.uint8) - 33].mean())
              for r in pysam.FastxFile(path)]
        out[name] = q_histogram(qs)
    return out


def best_amplicon(read, same_chrom):
    """Find the amplicon whose span best matches a read's aligned span.

    The match score is overlap / union of the two spans (Jaccard index): 1 for identical spans, 0 for disjoint
    ones. A read from a neighbouring, overlapping amplicon overlaps this one but scores low, because the union
    includes the parts that do not overlap.

    Args:
        read: a mapped pysam.AlignedSegment (its reference_start and reference_end are used).
        same_chrom: list of (name, start, end) of the amplicons on the read's chromosome, 0-based half-open.

    Returns:
        (name, score) of the best-matching amplicon, or (None, 0.0) when the read overlaps no amplicon. On a
        tie the first amplicon in the list wins. The caller decides whether the score is high enough
        (MIN_SPAN_MATCH).
    """
    rs, re_ = read.reference_start, read.reference_end
    best, score = None, 0.0
    for n, s, e in same_chrom:
        j = max(0, min(re_, e) - max(rs, s)) / (max(re_, e) - min(rs, s))
        if j > score:
            best, score = n, j
    return best, score


def counted(read):
    """Tell whether a read counts towards read noise.

    Args:
        read: a pysam.AlignedSegment.

    Returns:
        True for a primary alignment (neither secondary nor supplementary) with MAPQ >= MIN_MAPQ. Unmapped reads
        have MAPQ 0, so they are always False.
    """
    return not (read.is_secondary or read.is_supplementary or read.mapping_quality < MIN_MAPQ)


def amplicon_read_quality(bam, name, beds):
    """Histogram the read quality of the amplicon reads of one amplicon.

    A read is an amplicon read of this amplicon when best_amplicon() picks it and the match score is at least
    MIN_SPAN_MATCH. Amplicons overlap, so simply overlapping the amplicon would mix in reads from neighbouring
    amplicons. Only reads accepted by counted() are considered.

    Args:
        bam: open, indexed pysam.AlignmentFile.
        name: amplicon name, a key of beds.
        beds: {name: (chrom, start, end)} of every amplicon, 0-based half-open.

    Returns:
        The q_histogram() of the amplicon reads' mean qualities (computed from the BAM's base qualities, as in
        read_quality()), plus overlapping: the number of counted reads overlapping the amplicon at all.
    """
    chrom, start, end = beds[name]
    same_chrom = [(n, s, e) for n, (c, s, e) in beds.items() if c == chrom]
    qs, overlapping = [], 0
    for read in filter(counted, bam.fetch(chrom, start, end)):
        overlapping += 1
        best, score = best_amplicon(read, same_chrom)
        if best == name and score >= MIN_SPAN_MATCH:
            qs.append(float(-10 * np.log10(ERR_P[np.asarray(read.query_qualities)].mean())))
    return {**q_histogram(qs), "overlapping": overlapping}


def unmatched_reads(bam, beds):
    """Count the on-target reads and how many of them match no amplicon.

    A read is on target when counted() accepts it and it overlaps at least one amplicon; it matches no amplicon
    when its best_amplicon() score is below MIN_SPAN_MATCH (e.g. a chimeric or truncated product).

    Args:
        bam: open, indexed pysam.AlignmentFile.
        beds: {name: (chrom, start, end)} of every amplicon, 0-based half-open.

    Returns:
        {"on_target": int, "unmatched": int}.
    """
    on_target = unmatched = 0
    for chrom in sorted({c for c, _, _ in beds.values()}):
        same_chrom = [(n, s, e) for n, (c, s, e) in beds.items() if c == chrom]
        lo, hi = min(s for _, s, _ in same_chrom), max(e for _, _, e in same_chrom)
        for read in filter(counted, bam.fetch(chrom, lo, hi)):
            best, score = best_amplicon(read, same_chrom)
            if best is None:  # in a gap between amplicons
                continue
            on_target += 1
            unmatched += score < MIN_SPAN_MATCH
    return {"on_target": on_target, "unmatched": unmatched}


def overlapping_amplicons(name, beds, pools):
    """List the amplicons for the amplicons track of one amplicon's view.

    Args:
        name: the amplicon shown, a key of beds.
        beds: {name: (chrom, start, end)} of every amplicon, 0-based half-open.
        pools: {name: pool} from BED column 5; may be empty.

    Returns:
        One dict per amplicon (name, start, end as 1-based closed, pool or None): the amplicon shown first, then
        every other amplicon overlapping it by at least 1 bp, sorted by start.
    """
    chrom, start, end = beds[name]
    others = sorted((s, n) for n, (c, s, e) in beds.items() if c == chrom and n != name and s < end and e > start)
    return [{"name": n, "start": beds[n][1] + 1, "end": beds[n][2], "pool": pools.get(n)}
            for n in [name] + [n for _, n in others]]


def bam_summary(path):
    """Count the primary reads (mapped and unmapped) in a whole BAM.

    Note: the pipeline's primer-clipped BAMs miss the reads that clipping dropped, so this undercounts the
    dataset; the FASTQ is the reliable source for the dataset size (see read_quality()).

    Args:
        path: BAM path.

    Returns:
        {"path": path, "reads": number of primary reads}.
    """
    with pysam.AlignmentFile(path) as bam:
        n = sum(1 for r in bam.fetch(until_eof=True) if not (r.is_secondary or r.is_supplementary))
    return {"path": path, "reads": n}


def open_bam(path):
    """Open a BAM for region queries.

    Args:
        path: BAM path; its index (.bai) must exist.

    Returns:
        The open pysam.AlignmentFile.

    Raises:
        SystemExit: when the BAM has no index, with the samtools commands to make one.
    """
    bam =pysam.AlignmentFile(path)
    if not bam.has_index():
        raise SystemExit(f"{path} has no index; sort and index it first:\n"
                         f"  samtools sort -o sorted.bam {path} && samtools index sorted.bam")
    return bam


def main():
    """Build prototype_data.js (the global `DATA` read by prototype_layout.html) from the HD829 M20_4 files.

    Reads the BAM, VCF, FASTQ and BED paths set at the top of this file; command-line options add the known
    variants, an optional unfiltered BAM and a subset of amplicons. For each amplicon it stores the read noise of
    each dataset (count()), homopolymers, exons, overlapping amplicons, per-amplicon read quality and variants;
    for the sample it stores the dataset sizes and the read-quality histograms.
    """
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--known-variants", default=KNOWN_VARIANTS, help="CSV of known variants (default: %(default)s)")
    ap.add_argument("--unfiltered-bam", help="sorted, indexed BAM of the unfiltered reads, e.g. the reference workflow's "
                    "references/align_unfiltered_bam/*.unfiltered.aligned.sorted.bam (optional)")
    ap.add_argument("--amplicons", nargs="+", help="amplicons to include (default: every amplicon in the BED, in BED order)")
    args = ap.parse_args()
    known = known_variants(args.known_variants)
    bam = open_bam(BAM)
    ubam = open_bam(args.unfiltered_bam) if args.unfiltered_bam else None
    rows = [f for f in (line.split() for line in open(BED)) if len(f) >= 4]
    beds = {f[3]: (f[0], int(f[1]), int(f[2])) for f in rows}
    pools = {f[3]: f[4] for f in rows if len(f) >= 5}  # BED column 5 = pool
    datasets = {"filtered": {**bam_summary(BAM), **unmatched_reads(bam, beds)},
                "unfiltered": {**bam_summary(args.unfiltered_bam), **unmatched_reads(ubam, beds)} if ubam else None}
    data = {"sample": SAMPLE, "amplicons": [], "datasets": datasets, "read_quality": read_quality(),
            "params": {"min_mapq": MIN_MAPQ, "min_depth": MIN_DEPTH, "min_span_match": MIN_SPAN_MATCH}}
    unknown = set(args.amplicons or []) - set(beds)
    if unknown:
        raise SystemExit(f"not in {BED}: {', '.join(sorted(unknown))}")
    for name in args.amplicons or list(beds):
        chrom, start, end = beds[name]
        ref = ucsc_sequence(chrom, start, end)
        data["amplicons"].append({
            "name": name, "chrom": chrom, "start": start + 1, "end": end,
            "noise": count(bam, chrom, start, end, ref),
            "noise_unfiltered": count(ubam, chrom, start, end, ref) if ubam else None,
            "homopolymers": homopolymers(chrom, start, end),
            "exons": exons(chrom, start, end),
            "overlapping_amplicons": overlapping_amplicons(name, beds, pools),
            "read_quality": {"min_span_match": MIN_SPAN_MATCH, "filtered": amplicon_read_quality(bam, name, beds),
                             "unfiltered": amplicon_read_quality(ubam, name, beds) if ubam else None},
            "variants": variants(chrom, start, end, known),
        })
        print(name, "done", flush=True)
    with open(OUT, "w") as f:
        f.write("// PROTOTYPE data, generated by build_prototype_data.py\nconst DATA = ")
        json.dump(data, f, separators=(",", ":"))
        f.write(";\n")
    print("wrote", OUT, f"{os.path.getsize(OUT) / 1e6:.1f} MB")


if __name__ == "__main__":
    main()
