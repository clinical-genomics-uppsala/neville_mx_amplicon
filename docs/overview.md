# Pipeline overview

The pipeline has two entry points:

| Snakefile | Purpose |
|-----------|---------|
| `workflow/Snakefile` | Main analysis of one or more samples from a sequencing run. |
| `workflow/Snakefile_references` | Builds the artifact and background panels from a set of reference samples. See [Building reference panels](running_references.md). |

## Analysis steps

| Step | Tools | Page |
|------|-------|------|
| Basecalling, demultiplexing, adapter trimming | Dorado | [Basecalling](basecalling.md) |
| Quality and length filtering | Filtlong, Samtools | [Pre-alignment](prealignment.md) |
| Alignment, primer soft-clipping, split per target | Dorado (minimap2), Samtools | [Alignment](alignment.md) |
| SNV and indel calling, ensemble, annotation, filtering | ClairS-TO, DeepSomatic, VarDict, VT, bcbio-variation-recall, VEP, bcftools | [SNVs and indels](snv_indels.md) |
| Structural variant calling | Sniffles2 | [Structural variants](cnv_sv.md) |
| Quality control | pycoQC, Sequali, Mosdepth, Picard, MultiQC | [QC](qc.md) |
| Report | xlsxwriter | [Excel report](report.md) |

## Naming conventions

Intermediate files are named after the wildcards below.

| Wildcard | Meaning |
|----------|---------|
| `experiment` | Batch ID of the sequencing run (`batchid` in the config), e.g. `M21`. |
| `sample` | Sample ID from `samples.tsv`. |
| `type` | Sample type from `units.tsv`; `T` (tumour) for this analysis. |
| `target` | An amplicon (`amplicons` in the config) or a sub-region of an overlapping amplicon (`extra_regions`). |
| `pooln` | Pool number (keys of `pools` in the config). |
| `caller` | Variant caller: `clairs_to`, `deepsomatic`, `vardict` or `sniffles2`. |

Intermediate files are temporary by default; the final files are copied to `Results/` (see [Result files](result_files.md)).
Run Snakemake with `--notemp` to keep everything.

## Single-sample and multi-sample runs

The config key `multisample` switches the basecalling rules:

* `multisample: False`: one non-barcoded sample per flowcell (e.g. Flongle). POD5 files are basecalled directly into a sample BAM file.
* `multisample: True`: several barcoded samples on one flowcell (e.g. MinION). All POD5 files are basecalled together,
  then demultiplexed with a sample sheet, renamed per sample and trimmed.

Both paths produce the same file, `basecalling/dorado_duplex/{experiment}_{sample}_{type}_reads.ont_adapt_trim.bam`, and the
rest of the pipeline is identical.

!!! note
    In practice, the start scripts in `workflow/scripts/` often skip basecalling: they merge the BAM files already basecalled by
    MinKNOW (`bam_pass/`) and place them where the pipeline expects the basecalled BAM file. See [Running the pipeline](running.md).

## Rule graph

![rule graph](images/neville_mx_amplicon.png)

!!! note
    The rule graph image is generated from an earlier version of the pipeline, so some rule names differ from the current ones
    (for example `filtlong` is now `prealignment_ont_filtlong`). To regenerate it, edit the run paths in
    `workflow/scripts/rulegraph_multisample.sh` and run it from the repository root.
