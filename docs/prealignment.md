# Pre-alignment

## Read filtering with Filtlong

`prealignment_ont_filtlong` filters the basecalled reads with **[Filtlong](https://github.com/rrwick/Filtlong)** v0.2.1, without a reference:

```
filtlong --min_length 2000 --max_length 4000 --min_mean_q 30 --keep_percent 75
```

* **Length**: reads shorter than 2,000 bp are probably incomplete PCR products or reads truncated in the pore.
  Reads longer than 4,000 bp are probably chimeras of concatenated amplicons. Adjust both thresholds to the amplicon lengths of your design.
* **Quality**: reads with a mean quality below Q30 are discarded, and the worst reads are removed until 75% of the bases are kept.

Output: `prealignment/filtlong/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.fastq.gz`.

**Software settings**

| Option | Value | Description |
|--------|-------|-------------|
| `min_length` | 2000 | Minimum read length. |
| `max_length` | 4000 | Maximum read length. |

`--min_mean_q 30` and `--keep_percent 75` are set in the rule.

## Discarded reads

`prealignment_ont_fetch_filtered_reads` collects the reads removed by Filtlong into
`prealignment/filtlong/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.out.fastq.gz`.
Sequali reports their quality next to the kept reads (see [QC](qc.md#read-quality)).

`prealignment_ont_fetch_too_short_too_long_reads` splits the discarded reads into too-short and too-long reads
(`*.filtered.out.short.fastq.gz`, `*.filtered.out.long.fastq.gz`).
