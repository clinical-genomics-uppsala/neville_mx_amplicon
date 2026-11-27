# neville_mx_amplicon

**A Snakemake pipeline to basecall, filter, align and analyse Oxford Nanopore sequence data from pooled (multiplexed) amplicons.**

The pipeline is developed at Clinical Genomics Uppsala (CGU) to support rapid diagnostics of Acute Myeloid Leukemia (AML).
It is built on the [hydra-genetics](https://hydra-genetics.readthedocs.io/en/latest/) framework: generic steps are imported
from hydra-genetics modules, and the steps specific to ONT long reads and pooled amplicons are implemented in this repository.

Source code: [github.com/clinical-genomics-uppsala/neville_mx_amplicon](https://github.com/clinical-genomics-uppsala/neville_mx_amplicon)

## The assay

Multiplexed amplicon sequencing amplifies several targets together in a **multiplexed PCR**.
In the current design, 15 amplicons covering 5 genes (*IDH1*, *IDH2*, *NPM1*, *FLT3*, *TP53*) are amplified in 3 pools of 5 amplicons.
Overlapping amplicons are assigned to different pools to avoid truncated PCR products.
Several barcoded samples can also be sequenced together on the same flowcell.
The libraries are sequenced on a MinION device with ONT sequencing kit 14 (R10.4.1 flowcells).

| Pool | Amplicon | Length (bp) | Chromosome |
|------|----------|-------------|------------|
| 1 | TP53_A | 3000 | chr17 |
| 1 | TP53_D | 3000 | chr17 |
| 1 | TP53_F | 3000 | chr17 |
| 1 | TP53_I | 3000 | chr17 |
| 1 | FLT3_TKD | 3000 | chr13 |
| 2 | IDH1 | 2000 | chr2 |
| 2 | TP53_B | 3000 | chr17 |
| 2 | TP53_G | 3000 | chr17 |
| 2 | TP53_J | 3000 | chr17 |
| 2 | NPM1 | 2000 | chr5 |
| 3 | IDH2 | 2000 | chr15 |
| 3 | TP53_C | 3000 | chr17 |
| 3 | TP53_E | 3000 | chr17 |
| 3 | TP53_H | 3000 | chr17 |
| 3 | FLT3_ITD | 3000 | chr13 |

The whole *TP53* gene is tiled by 10 overlapping amplicons (A to J). Amplicons TP53_D and TP53_J overlap completely, which needs
special handling when estimating per-amplicon read counts (see [Alignment](alignment.md#tp53-j-amplicon)).

## What the pipeline does

1. **Basecalling** of POD5 files with Dorado, with demultiplexing and adapter trimming for barcoded runs.
2. **Read filtering** on quality and length with Filtlong.
3. **Alignment** with Dorado (minimap2), then **soft-clipping of primers** with `samtools ampliconclip`.
4. **Small variant calling** with ClairS-TO, DeepSomatic and VarDict, merged into one ensemble VCF.
5. **Structural variant calling** with Sniffles2.
6. **Annotation** with VEP, **artifact annotation** against a panel of reference samples, and soft filtering.
7. **QC**: read quality, coverage per amplicon, exon and pool, on-target rate, all summarised in MultiQC.
8. **Excel report** with filtered variants and coverage tables per sample.

See the [pipeline overview](overview.md) for the full description and the rule graph.

## Where to start

* New user: read [Pipeline setup and configuration](setup.md), then [Running the pipeline](running.md).
* To try the pipeline on a laptop, run the [small integration test](integration_test.md).
* To find a specific output file, see [Result files](result_files.md).
* To build artifact and background panels from normal samples, see [Building reference panels](running_references.md).

## Hydra-genetics

[Hydra-genetics](https://github.com/hydra-genetics) is a community that aims to make
[Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) pipeline development easier, faster, more structured and of higher quality.
It provides tested [Snakemake modules](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules)
that can be combined into complete analyses. This pipeline imports the following modules:

| Module | Version | Used for |
|--------|---------|----------|
| [snv_indels](https://hydra-genetics-snv-indels.readthedocs.io/en/latest/) | v1.1.0 | VarDict, VT decompose/normalize, ensemble calling, bgzip/tabix |
| [cnv_sv](https://hydra-genetics-cnv-sv.readthedocs.io/en/latest/) | v0.9.0 | Sniffles2 |
| [annotation](https://hydra-genetics-annotation.readthedocs.io/en/latest/) | v1.0.0 | VEP |
| [filtering](https://hydra-genetics-filtering.readthedocs.io/en/latest/) | v1.1.0 | Soft filtering of VCF files |
| [qc](https://hydra-genetics-qc.readthedocs.io/en/latest/) | v0.5.0 | Mosdepth, Picard, MultiQC |
| [references](https://hydra-genetics-references.readthedocs.io/en/latest/) | v1.0.1 | Background panels (reference workflow only) |

Module versions are set under `modules:` in `config/config.yaml`.
