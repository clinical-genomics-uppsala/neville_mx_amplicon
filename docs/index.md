# Welcome to neville_mx_amplicon

`neville_mx_amplicon` is a Snakemake-based bioinformatics pipeline designed for the filtering, alignment, variant calling, and quality control of Oxford Nanopore Technologies (ONT) sequencing data generated from pooled multiplexed amplicons. 

Developed at Clinical Genomics Uppsala (CGU), this workflow is tailored to support rapid diagnostic settings, specifically targeting mutations relevant in Acute Myeloid Leukemia (AML) diagnostics.

## Features

- **Basecalling & Demultiplexing**: Uses *Dorado* to process raw POD5 signal data and perform barcode demultiplexing.
- **Read Quality Filtering**: Employs *Filtlong* to filter reads based on length and quality threshold constraints.
- **Alignment & Processing**: Utilizes *Dorado* and *Samtools* to align reads and perform primer soft-clipping.
- **Variant Calling**: Incorporates multiple long-read somatic and structural variant callers:
    - *ClairS-TO*
    - *DeepSomatic*
    - *VarDict*
    - *Sniffles2*
- **Annotation & Post-processing**: Normalizes variants using *VT* and annotates using *VEP*.
- **Quality Control**: Generates QC reports using *Sequali*, *Picard*, *Mosdepth*, and compiles them into interactive *MultiQC* dashboards and Excel reports.

---

## Hydra-Genetics Community

This pipeline is built on the [hydra-genetics](https://github.com/hydra-genetics) framework, a community-driven initiative with the goal of making Snakemake pipeline development easier, faster, more structured, and of higher quality.

All modules are subjected to extensive testing to make sure that new releases do not unexpectedly break existing pipelines or deviate from development guidelines.
