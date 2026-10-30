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
AMPLICONS = ["TP53_D2", "IDH1_1"]
MIN_MAPQ, MIN_DEPTH, MIN_SPAN_MATCH = 20, 100, 0.5
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "prototype_data.js")

LUT = np.full(256, 4, dtype=np.int64)
for i, b in enumerate(b"ACGT"):
    LUT[b] = i


def ucsc_sequence(chrom, start, end):
    url = f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom={chrom};start={start};end={end}"
    with urllib.request.urlopen(url) as r:
        return json.load(r)["dna"].upper()


def exons(chrom, start, end):
    url = f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=ncbiRefSeqSelect;chrom={chrom};start={start};end={end}"
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
    """Known variants of the reference standard. Required columns: chrom, pos, ref, alt (GRCh38, VCF-style alleles).
    Optional: expected_af, gene, hgvsc, hgvsp, note. Comma, semicolon or tab delimited."""
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
        name = " ".join(x for x in (r.get("gene"), r.get("hgvsp") or r.get("hgvsc")) if x and x != ".")
        r["label"] = name or f"{r['chrom']}:{r['pos']} {r['ref']}>{r['alt']}"
    return [r for r in rows if r["chrom"] and r["pos"]]


def count(bam, chrom, start, end, ref):
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
    r3 = lambda a: [round(float(v), 5) for v in a]  # noqa: E731
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
    hp = []
    for line in open(HOMOPOLYMERS):
        c, s, e, base, ln = line.split()[:5]
        if c == chrom and int(s) >= start and int(e) <= end:
            hp.append({"start": int(s) + 1, "end": int(e), "base": base, "length": int(ln)})
    return hp


ERR_P = 10 ** (-np.arange(94) / 10)
Q_BINS = np.arange(0, 50.5, 0.5)


def q_histogram(qs):
    h, _ = np.histogram(qs, bins=Q_BINS)
    return {"counts": h.tolist(), "n": len(qs), "median": round(float(np.median(qs)), 2) if qs else None}


def read_quality():
    out = {"bins": Q_BINS[:-1].tolist()}
    for name, path in FASTQ.items():
        qs = [-10 * np.log10(ERR_P[np.frombuffer(r.quality.encode(), dtype=np.uint8) - 33].mean())
              for r in pysam.FastxFile(path)]
        out[name] = q_histogram(qs)
    return out


def best_amplicon(read, same_chrom):
    """The amplicon whose span best matches the read's aligned span (overlap / union), and that score."""
    rs, re_ = read.reference_start, read.reference_end
    best, score = None, 0.0
    for n, s, e in same_chrom:
        j = max(0, min(re_, e) - max(rs, s)) / (max(re_, e) - min(rs, s))
        if j > score:
            best, score = n, j
    return best, score


def counted(read):
    return not (read.is_secondary or read.is_supplementary or read.mapping_quality < MIN_MAPQ)


def amplicon_read_quality(bam, name, beds):
    """Read quality of the amplicon reads of one amplicon: a read belongs to the amplicon whose span best matches
    its aligned span (overlap / union), if that ratio is at least MIN_SPAN_MATCH. Amplicons overlap, so plain
    overlap with the amplicon would mix in reads from neighbouring amplicons."""
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
    """Reads overlapping at least one amplicon (primary, MAPQ >= MIN_MAPQ) and how many of them match no amplicon,
    i.e. their best overlap / union is below MIN_SPAN_MATCH."""
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
    """Every amplicon of the BED overlapping this one by at least 1 bp, this one first (1-based closed spans)."""
    chrom, start, end = beds[name]
    others = sorted((s, n) for n, (c, s, e) in beds.items() if c == chrom and n != name and s < end and e > start)
    return [{"name": n, "start": beds[n][1] + 1, "end": beds[n][2], "pool": pools.get(n)}
            for n in [name] + [n for _, n in others]]


def bam_summary(path):
    """Primary reads in the whole BAM (mapped and unmapped), i.e. the size of the dataset."""
    with pysam.AlignmentFile(path) as bam:
        n = sum(1 for r in bam.fetch(until_eof=True) if not (r.is_secondary or r.is_supplementary))
    return {"path": path, "reads": n}


def open_bam(path):
    bam = pysam.AlignmentFile(path)
    if not bam.has_index():
        raise SystemExit(f"{path} has no index; sort and index it first:\n"
                         f"  samtools sort -o sorted.bam {path} && samtools index sorted.bam")
    return bam


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--known-variants", default=KNOWN_VARIANTS, help="CSV of known variants (default: %(default)s)")
    ap.add_argument("--unfiltered-bam", help="sorted, indexed BAM of the unfiltered reads, e.g. the reference workflow's "
                    "references/align_unfiltered_bam/*.unfiltered.aligned.sorted.bam (optional)")
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
    for name in AMPLICONS:
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
