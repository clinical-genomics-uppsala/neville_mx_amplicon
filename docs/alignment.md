# Alignment

## Alignment with Dorado

`alignment_ont_dorado_align` aligns the filtered reads to the reference genome with `dorado aligner`, which uses
**[minimap2](https://github.com/lh3/minimap2)** with the default ONT preset (`lr:hq`).
The BAM file is then sorted (`alignment_ont_bam_sort`) and indexed (`alignment_ont_bam_index`) with **Samtools** v1.21.

**Software settings**

| Option | Value | Description |
|--------|-------|-------------|
| `extra` | | Extra options for `dorado aligner`, e.g. `--mm2-opts "-x sr"`. |

## Primer soft-clipping

`alignment_ont_bam_softclip` soft-clips the primers with
**[samtools ampliconclip](https://www.htslib.org/doc/samtools-ampliconclip.html)** and the primer BED file `amplicons-primers.bed`:

```
samtools ampliconclip --soft-clip --both-ends --clipped -b amplicons-primers.bed
```

Primer bases come from the oligonucleotides, not from the sample, so they must not be used for variant calling.
Soft-clipping keeps them in the reads but excludes them from pileups.
The summary of clipped reads (`results/qc/samtools_ampliconclip/`) is shown in MultiQC.

Output: `alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam`.
This is the BAM file used for variant calling and QC, and delivered in `Results/Data/`.

## Split per target

`alignment_ont_split_bam_by_target` writes one BAM file per target (each entry in `amplicons` and `extra_regions`) with
`samtools view -L {target}.bed`, plus the list of read names in that BAM file.
These files are used for per-target coverage and alignment metrics.

### TP53 J amplicon

TP53_D2 and TP53_J3 overlap completely, so there is no region covered only by J3.
`alignment_ont_split_bam_target_j3` builds a `TP53_J3_only` BAM file by taking the reads in the `TP53_D2+J3` region and
removing the reads found in the `TP53_D2_only` region.
In the same way, the read count of J3 in each pool is the count for `TP53_D2+J3` minus the count for `TP53_D2_only`
(see [Yield per pool](qc.md#yield-per-amplicon-and-per-pool)).
