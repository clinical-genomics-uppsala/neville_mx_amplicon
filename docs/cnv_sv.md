# Structural variants

Structural variants (SVs), such as the *FLT3* internal tandem duplications and the *NPM1* insertions, are called with
**[Sniffles2](https://github.com/fritzsedlazeck/Sniffles)** v2.6.1 from the hydra-genetics
[cnv_sv](https://hydra-genetics-cnv-sv.readthedocs.io/en/latest/) module.

## Pipeline output files

* `Results/Data/{experiment}_{sample}_{type}.ont_reads.cnv-sv.vcf.gz`
* SVs sheet of the [Excel report](report.md)

## Sniffles2 calling

`cnv_sv_sniffles2_call` runs on the soft-clipped BAM file.

| Option | Value | Description |
|--------|-------|-------------|
| `extra` | `--minsvlen 3` | Report SVs from 3 bp, so that short insertions such as *NPM1* type A (4 bp) are called. |

## Region filtering

Sniffles2 cannot be restricted to regions, so `prefiltering_bcftools_include_region` keeps only the calls in
`amplicons-sniffles2.bed` with `bcftools filter -R` (set in `prefiltering_bcftools_include_region.filter.sniffles2.extra`).
The filtered VCF is annotated with VEP as for the small variants, but there is no artifact annotation or soft filtering for SVs.

## Allele frequency

Sniffles2 reports `VAF` in the INFO field, computed as
[`SUPPORT / COVERAGE`](https://github.com/fritzsedlazeck/Sniffles/blob/7113f927c3f9ca6329c436cd25cec369a100da31/src/sniffles/genotyping.py#L90),
where `SUPPORT` is the number of reads supporting the SV and `COVERAGE` depends on the SV type:

* **insertions**: coverage at the centre of the insertion;
* **duplications**: mean of the coverage at the start and the end, plus 75% of the support;
* **inversions**: mean of the coverage upstream and downstream, plus 50% of the support;
* **deletions**: mean of the coverage at the start, centre and end, each increased by the support from supplementary alignments.

## Other callers

[ScanITD](https://github.com/ylab-hi/ScanITD) (`cnv_sv_scanitd`) and SAVANA are kept in the code but disabled.
