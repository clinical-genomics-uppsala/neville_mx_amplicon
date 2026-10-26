# HD829 known variants in GRCh38

Research note for issue #36, part of map #35 (noise visualisation tool for ONT amplicon data).
This branch is a throwaway record and will not be merged.

## Answer

- HD829 is Horizon Discovery's **Mimix Myeloid Cancer Panel, gDNA Reference Standard**, also called the "Myeloid gDNA reference standard HD829". It contains **22 variants** in 19 genes: 15 at 5 % expected allele frequency (AF) and 7 at 10 to 70 %.
- **21 of the 22** variants resolve unambiguously to GRCh38 chrom/pos/ref/alt (16 SNVs, 5 small indels). The one that does not is **FLT3 ITD300**, an engineered insertion of about 300 bp. Horizon gives only its position (13:28033909), not the inserted sequence.
- **6 variants** fall inside the 15 amplicons of `HD829/M20_4/Bed/amplicons.bed`. They are FLT3 D835Y, FLT3 ITD300, IDH1 R132C, IDH2 R172K, NPM1 W288Cfs\*12 and TP53 S241F. In the M20_4 ONT run, only TP53 S241F was called (FILTER=PASS, AF=0.055).
- Proposed known-variants TSV: `chrom pos ref alt expected_af gene hgvsc hgvsp note`. See [Proposed TSV layout](#proposed-tsv-layout). The filled-in file is [`hd829_known_variants.tsv`](hd829_known_variants.tsv).

## Sources

| # | Source | What it gives | Accessed |
|---|---|---|---|
| S1 | Horizon product page, <https://horizondiscovery.com/en/reference-standards/products/myeloid-gdna-reference-standard> (page title "Mimix™ Myeloid Cancer Panel, gDNA Reference Standard") | Variant table: gene, Variant (AA), CDS mutation, GRCh38 co-ordinates, expected AF, NGS AF and read depth for one batch. Also the GATA1 AF update note. | 2026-10-26 |
| S2 | Horizon product data PDF "Shelf life extension for: Myeloid gDNA reference standard HD829" (Doc ID 03-02-013I), <https://horizondiscovery.com/-/media/Files/Horizon/resources/Product-data/hd829-shelf-life-extension.pdf?sc_lang=en>, linked from S1 | Ties product code **HD829** to the "Myeloid gDNA reference standard" sold on S1 | 2026-10-26 |
| S3 | Ensembl REST `POST /sequence/region/human` and `POST /vep/human/region?hgvs=1;mane=1` (GRCh38) | GRCh38 reference base at each S1 position. For SNVs: the alt allele whose HGVSc on the gene's transcripts equals S1's "CDS mutation", plus the MANE Select transcript. | 2026-10-26 |
| S4 | NCBI Variation Services: `/variation/v0/hgvs/{NM}:{c.}/contextuals`, `/spdi/{spdi}/all_equivalent_contextual`, `/spdi/{spdi}/vcf_fields` | For indels: maps the MANE RefSeq c. change to a GRCh38 (NC_0000xx) SPDI, then to left-normalised VCF fields | 2026-10-26 |
| S5 | Local `/home/camille/ampliconthemato/HD829/HD829_47035.hg38.genes.dbsnp.vcf` | Independent check. All 21 resolved variants appear there with identical pos/ref/alt, and its DP values equal S1's "NGS read depth" column exactly. This file therefore looks like Horizon's exome data for batch 47035, but how it was obtained is not documented locally. | 2026-10-26 |

The Ensembl `variant_recoder` endpoint with gene-symbol HGVS (for example `TP53:c.722C>T`) timed out or returned 500 for nearly every query on the access date. That is why S3/S4 were used instead. The Horizon literature-review PDF linked from S1 returned HTML rather than a PDF, so it was not read.

**Is S1 the HD829 product?** Yes. `HD829/README.md` describes HD829 as a control sample used routinely in the tM pipeline, with "most [variants] at ca. 5 % VAF". S1 states that 15 of the 22 variants are at 5 %. S1 links the HD829 shelf-life document (S2), which names the product "Myeloid gDNA reference standard HD829". The order code itself is not printed in S1's static HTML. S5's depths match S1's table exactly.

## Method

1. Copied the variant table from S1 (raw HTML, not a rendered summary).
2. **SNVs (16):** fetched the GRCh38 reference base at S1's position (S3). Ran VEP on all three possible alt alleles and kept the one whose `hgvsc` equals S1's CDS change. All 16 matched exactly one alt, always on the gene's MANE Select transcript except for ABL1 (below).
3. **Indels (5):** expressed S1's c. change on the gene's MANE Select RefSeq transcript (taken from step 2) and mapped it to GRCh38 through NCBI SPDI (S4). Converted to left-normalised VCF fields (S4). In all 5 cases the VCF POS equals S1's coordinate, so S1 reports the VCF anchor base. GATA2 was also confirmed with Ensembl `variant_recoder` on `NM_032638.5:c.599del` → `3-128485998-AC-A`.
4. Cross-checked all 21 against S5: 21/21 identical.
5. Amplicon overlap: a variant counts as inside an amplicon when `bed_start < POS` and `POS + len(REF) - 1 <= bed_end` (BED is 0-based, half-open).
6. ONT cross-check (sanity only): looked for records within ±10 bp in `HD829/M20_4/Data/HD829_T.ont_reads.snv-indels.vcf` (ClairS-TO/DeepVariant output, annotated).

## Variants

GRCh38, VCF-style (1-based, left-normalised, anchor base for indels). "In amplicon" refers to `M20_4/Bed/amplicons.bed`.

| Gene | Horizon AA | Transcript : c. (MANE Select unless noted) | chrom | pos | ref | alt | Expected AF | In amplicon | M20_4 ONT VCF |
|---|---|---|---|---|---|---|---|---|---|
| ABL1 | T315I | NM_007313.3 / ENST00000372348.9 : c.1001C>T (p.T334I, non-MANE) | chr9 | 130872896 | C | T | 5 % | – | |
| ASXL1 | G646Wfs\*12 | NM_015338.6 : c.1934dup | chr20 | 32434638 | A | AG | 40 % | – | |
| ASXL1 | W796C | NM_015338.6 : c.2388G>T | chr20 | 32435100 | G | T | 5 % | – | |
| BCOR | Q1208Tfs\*8 | NM_001123385.2 : c.3621dup | chrX | 40063833 | G | GT | 70 % | – | |
| CBL | S403F | NM_005188.4 : c.1208C>T | chr11 | 119278278 | C | T | 5 % | – | |
| DNMT3A | R882C | NM_022552.5 : c.2644C>T | chr2 | 25234374 | G | A | 5 % | – | |
| EZH2 | R418Q | NM_004456.5 : c.1253G>A | chr7 | 148817379 | C | T | 5 % | – | |
| FLT3 | D835Y | NM_004119.3 : c.2503G>T | chr13 | 28018505 | C | A | 5 % | FLT3_TKD_3kb_1 | not called |
| FLT3 | ITD300 | – (no CDS given) | chr13 | 28033909 | C | **unresolved** (~300 bp insertion) | 5 % | FLT3_ITD_3kb_1 | not called |
| GATA1 | Q119\* | NM_002049.4 : c.355C>T | chrX | 48791978 | C | T | 11 % (\*) | – | |
| GATA2 | G200Vfs\*18 | NM_032638.5 : c.599del | chr3 | 128485998 | AC | A | 35 % | – | |
| IDH1 | R132C | NM_005896.4 : c.394C>T | chr2 | 208248389 | G | A | 5 % | IDH1_1 | not called |
| IDH2 | R172K | NM_002168.4 : c.515G>A | chr15 | 90088606 | C | T | 5 % | IDH2_1 | not called |
| JAK2 | F537-K539>L | NM_004972.4 : c.1611_1616del | chr9 | 5070021 | TTCACAA | T | 5 % | – | |
| JAK2 | V617F | NM_004972.4 : c.1849G>T | chr9 | 5073770 | G | T | 5 % | – | |
| KRAS | G13D | NM_004985.5 : c.38G>A | chr12 | 25245347 | C | T | 40 % | – | |
| NPM1 | W288Cfs\*12 | NM_002520.7 : c.860_863dup | chr5 | 171410539 | C | CTCTG | 5 % | NPM1_5 | not called |
| NRAS | Q61L | NM_002524.5 : c.182A>T | chr1 | 114713908 | T | A | 10 % | – | |
| RUNX1 | M267I | NM_001754.5 : c.801G>A | chr21 | 34834414 | C | T | 35 % | – | |
| SF3B1 | G740E | NM_012433.4 : c.2219G>A | chr2 | 197401989 | C | T | 5 % | – | |
| TET2 | R1261H | NM_001127208.3 : c.3782G>A | chr4 | 105243757 | G | A | 5 % | – | |
| TP53 | S241F | NM_000546.6 : c.722C>T | chr17 | 7674241 | G | A | 5 % | TP53_D2, TP53_J3, TP53_E1 | PASS, AF=0.055, DP=6423 |

(\*) S1 says the GATA1 Q119\* expected AF was raised from 10 % to 11 %. Batch 67192 and all batches after 67533 report 11 % on the Certificate of Analysis. Which batch was sequenced as HD829 in M20_4 is not recorded in this repo. The batch-specific Certificate of Analysis is the authority for exact AFs.

"Not called" means no record within ±10 bp in the M20_4 VCF. That VCF holds calls only on chr2, chr5, chr13, chr15 and chr17, which matches the amplicon panel. The missing calls are consistent with `HD829/README.md`, which says most expected variants were not called.

## Caveats and unresolved items

1. **FLT3 ITD300 is unresolved.** Horizon describes it as an engineered, novel 300 bp insertion. S1 gives only `13:28033909` and no inserted sequence ("N/A" CDS), and S1's own whole-exome table did not call it. The TSV carries it as `ref=C` (GRCh38 base at that position) and `alt=<INS>` (VCF symbolic allele) with a note. A tool should treat it as "look for a ~300 bp insertion near this position", not as an exact allele. The exact sequence would have to come from Horizon or from de novo assembly of the FLT3_ITD amplicon reads.
2. **ABL1 naming.** S1 lists "T315I" (the conventional isoform 1a name; NM_005157 c.944C>T) together with "c.1001C>T", which is isoform 1b numbering (NM_007313.3, ENST00000372348; protein p.T334I). Both describe the same genomic change, chr9:130872896 C>T, so the coordinate is unambiguous. Only the label differs.
3. **Expected AF is per batch.** S1's "Expected allelic frequency" is the product specification. The table's NGS AF and depth come from one example batch. The Certificate of Analysis for the lot actually used is not available here.
4. Transcript choice does not affect any coordinate. For every SNV, all transcripts that matched gave the same genomic allele, and S5 agrees for all 21.

## Proposed TSV layout

The layout is generic enough for any reference standard (Horizon, SeraCare, in-house mixes): one tab-separated file per standard, one header row, one row per expected variant.

| Column | Required | Meaning |
|---|---|---|
| `chrom` | yes | Contig name exactly as in the reference FASTA/BAM (for example `chr13`) |
| `pos` | yes | 1-based VCF POS (anchor base for indels) |
| `ref` | yes | VCF REF, left-normalised against the same reference |
| `alt` | yes | VCF ALT; a symbolic allele (`<INS>`, `<DEL>`, `<DUP>`) is allowed when the exact sequence is unknown |
| `expected_af` | yes | Expected variant allele fraction, 0–1 (not percent) |
| `gene` | no | HGNC symbol, for labels |
| `hgvsc` | no | `transcript:c.` change, transcript version included (for example `NM_000546.6:c.722C>T`) |
| `hgvsp` | no | Protein change used for display (vendor's label may differ, see `note`) |
| `note` | no | Free text: provenance, batch caveats, unresolved status |

Why this shape:

- `chrom/pos/ref/alt` follow VCF conventions, so rows can be matched directly against caller output and BAM pileups after normalisation (for example `bcftools norm`). They also turn into a VCF with a trivial header.
- AF is a fraction, matching FORMAT/AF in the pipeline's VCFs.
- Only the first five columns are needed to plot expected vs observed. The rest are labels.
- The genome build is not a column. It belongs in the file name or a `#` comment line (for example `# build=GRCh38`), so it is not repeated on every row.
