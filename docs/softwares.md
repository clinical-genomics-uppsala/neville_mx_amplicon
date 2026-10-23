# Pipeline specific rules

Rules defined in `workflow/rules/`. Rules imported from hydra-genetics modules are documented in the
[hydra-genetics module documentation](https://hydra-genetics.readthedocs.io/en/latest/).

# Basecalling

Rules in `workflow/rules/basecalling.smk`.

## [basecalling_dorado_basecaller](https://software-docs.nanoporetech.com/dorado/latest/basecaller/simplex/)
Simplex basecalling of POD5 files with Dorado, with adapter trimming (single-sample runs).

### :snake: Rule

This rule is defined inside an `if config['multisample']` block; see its source in [workflow/rules/basecalling.smk](https://github.com/clinical-genomics-uppsala/neville_mx_amplicon/blob/develop/workflow/rules/basecalling.smk).

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__basecalling_dorado_basecaller#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__basecalling_dorado_basecaller#


## [basecalling_dorado_duplex](https://software-docs.nanoporetech.com/dorado/latest/basecaller/duplex/)
Duplex basecalling of POD5 files with Dorado, piped to `dorado trim` (single-sample runs).

### :snake: Rule

This rule is defined inside an `if config['multisample']` block; see its source in [workflow/rules/basecalling.smk](https://github.com/clinical-genomics-uppsala/neville_mx_amplicon/blob/develop/workflow/rules/basecalling.smk).

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__basecalling_dorado_duplex#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__basecalling_dorado_duplex#


## [basecalling_dorado_duplex_multisamples](https://software-docs.nanoporetech.com/dorado/latest/basecaller/duplex/)
Duplex basecalling of all POD5 files of a multiplexed run, without trimming so that barcodes are kept.

### :snake: Rule

This rule is defined inside an `if config['multisample']` block; see its source in [workflow/rules/basecalling.smk](https://github.com/clinical-genomics-uppsala/neville_mx_amplicon/blob/develop/workflow/rules/basecalling.smk).

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__basecalling_dorado_duplex_multisamples#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__basecalling_dorado_duplex_multisamples#


## [basecalling_dorado_demux](https://software-docs.nanoporetech.com/dorado/latest/barcoding/barcoding/)
Demultiplexing of the basecalled reads with `dorado demux` and the sample sheet.

### :snake: Rule

This rule is defined inside an `if config['multisample']` block; see its source in [workflow/rules/basecalling.smk](https://github.com/clinical-genomics-uppsala/neville_mx_amplicon/blob/develop/workflow/rules/basecalling.smk).

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__basecalling_dorado_demux#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__basecalling_dorado_demux#


## basecalling_rename_demux_bam
Copy each demultiplexed BAM file to a file named after the sample.

### :snake: Rule

This rule is defined inside an `if config['multisample']` block; see its source in [workflow/rules/basecalling.smk](https://github.com/clinical-genomics-uppsala/neville_mx_amplicon/blob/develop/workflow/rules/basecalling.smk).

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__basecalling_rename_demux_bam#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__basecalling_rename_demux_bam#


## [basecalling_dorado_trim](https://software-docs.nanoporetech.com/dorado/latest/basecaller/read_trimming/)
Adapter and barcode trimming of demultiplexed reads with `dorado trim`.

### :snake: Rule

This rule is defined inside an `if config['multisample']` block; see its source in [workflow/rules/basecalling.smk](https://github.com/clinical-genomics-uppsala/neville_mx_amplicon/blob/develop/workflow/rules/basecalling.smk).

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__basecalling_dorado_trim#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__basecalling_dorado_trim#


## [basecalling_bam_to_fastq](https://www.htslib.org/doc/samtools-fasta.html)
Convert the unaligned BAM file to FASTQ with `samtools fastq`.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__basecalling__basecalling_bam_to_fastq#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__basecalling__basecalling_bam_to_fastq#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__basecalling_bam_to_fastq#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__basecalling_bam_to_fastq#


## basecalling_compress_fastq
Compress the FASTQ file with gzip.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__basecalling__basecalling_compress_fastq#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__basecalling__basecalling_compress_fastq#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__basecalling_compress_fastq#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__basecalling_compress_fastq#


# Pre-alignment

Rules in `workflow/rules/prealignment_ont.smk`.

## [prealignment_ont_filtlong](https://github.com/rrwick/Filtlong)
Filter reads on mean quality and length with Filtlong.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__prealignment_ont__prealignment_ont_filtlong#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__prealignment_ont_filtlong#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__prealignment_ont_filtlong#


## [prealignment_ont_fetch_filtered_reads](https://www.htslib.org/doc/samtools-view.html)
Collect the reads removed by Filtlong.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__prealignment_ont__prealignment_ont_fetch_filtered_reads#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__prealignment_ont__prealignment_ont_fetch_filtered_reads#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__prealignment_ont_fetch_filtered_reads#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__prealignment_ont_fetch_filtered_reads#


## prealignment_ont_fetch_too_short_too_long_reads
Split the removed reads into too-short and too-long reads.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__prealignment_ont__prealignment_ont_fetch_too_short_too_long_reads#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__prealignment_ont__prealignment_ont_fetch_too_short_too_long_reads#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__prealignment_ont_fetch_too_short_too_long_reads#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__prealignment_ont_fetch_too_short_too_long_reads#


# Alignment

Rules in `workflow/rules/alignment_ont.smk`.

## [alignment_ont_dorado_align](https://software-docs.nanoporetech.com/dorado/latest/basecaller/alignment/)
Align the filtered reads with `dorado aligner` (minimap2).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__alignment_ont__alignment_ont_dorado_align#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__alignment_ont__alignment_ont_dorado_align#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__alignment_ont_dorado_align#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__alignment_ont_dorado_align#


## [alignment_ont_bam_sort](https://www.htslib.org/doc/samtools-sort.html)
Sort the aligned reads.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__alignment_ont__alignment_ont_bam_sort#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__alignment_ont_bam_sort#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__alignment_ont_bam_sort#


## [alignment_ont_bam_index](https://www.htslib.org/doc/samtools-index.html)
Index the sorted BAM file.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__alignment_ont__alignment_ont_bam_index#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__alignment_ont_bam_index#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__alignment_ont_bam_index#


## [alignment_ont_bam_softclip](https://www.htslib.org/doc/samtools-ampliconclip.html)
Soft-clip the primers with `samtools ampliconclip`.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__alignment_ont__alignment_ont_bam_softclip#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__alignment_ont__alignment_ont_bam_softclip#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__alignment_ont_bam_softclip#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__alignment_ont_bam_softclip#


## [alignment_ont_split_bam_by_target](https://www.htslib.org/doc/samtools-view.html)
Write one BAM file and one list of read names per target.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__alignment_ont__alignment_ont_split_bam_by_target#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__alignment_ont__alignment_ont_split_bam_by_target#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__alignment_ont_split_bam_by_target#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__alignment_ont_split_bam_by_target#


## [alignment_ont_split_bam_target_j3](https://www.htslib.org/doc/samtools-view.html)
Build the BAM file for TP53_J3 from the D2+J3 region minus the reads of D2 only.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__alignment_ont__alignment_ont_split_bam_target_j3#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__alignment_ont__alignment_ont_split_bam_target_j3#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__alignment_ont_split_bam_target_j3#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__alignment_ont_split_bam_target_j3#


# Variant calling

Rules in `workflow/rules/snv_indels_svs.smk`.

## [snv_indels_svs_clairs_to](https://github.com/HKU-BAL/ClairS-TO)
Tumour-only small variant calling with ClairS-TO.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__snv_indels_svs__snv_indels_svs_clairs_to#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__snv_indels_svs__snv_indels_svs_clairs_to#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__snv_indels_svs_clairs_to#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__snv_indels_svs_clairs_to#


## [snv_indels_svs_clairs_to_concat](https://samtools.github.io/bcftools/bcftools.html#concat)
Concatenate and sort the SNV and indel VCF files from ClairS-TO.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__snv_indels_svs__snv_indels_svs_clairs_to_concat#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__snv_indels_svs__snv_indels_svs_clairs_to_concat#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__snv_indels_svs_clairs_to_concat#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__snv_indels_svs_clairs_to_concat#


## [snv_indels_svs_deepsomatic](https://github.com/google/deepsomatic)
Tumour-only small variant calling with DeepSomatic (`ONT_TUMOR_ONLY` model).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__snv_indels_svs__snv_indels_svs_deepsomatic#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__snv_indels_svs__snv_indels_svs_deepsomatic#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__snv_indels_svs_deepsomatic#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__snv_indels_svs_deepsomatic#


## [snv_indels_svs_savana](https://github.com/cortes-ciriano-lab/savana)
Somatic SV calling with SAVANA. Not used: SAVANA does not handle amplicon data yet.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__snv_indels_svs__snv_indels_svs_savana#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__snv_indels_svs__snv_indels_svs_savana#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__snv_indels_svs_savana#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__snv_indels_svs_savana#


# Pre-filtering and annotation

Rules in `workflow/rules/prefiltering.smk`.

## [prefiltering_bcftools_view](https://samtools.github.io/bcftools/bcftools.html#view)
Remove each caller's false-positive calls with caller-specific `bcftools view` expressions.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__prefiltering__prefiltering_bcftools_view#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__prefiltering__prefiltering_bcftools_view#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__prefiltering_bcftools_view#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__prefiltering_bcftools_view#


## [prefiltering_bcftools_include_region](https://samtools.github.io/bcftools/bcftools.html#filter)
Keep only the Sniffles2 calls in the target regions.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__prefiltering__prefiltering_bcftools_include_region#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__prefiltering__prefiltering_bcftools_include_region#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__prefiltering_bcftools_include_region#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__prefiltering_bcftools_include_region#


## [prefiltering_rename_vaf_to_af](https://samtools.github.io/bcftools/bcftools.html#annotate)
Rename `FORMAT/VAF` to `FORMAT/AF`.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__prefiltering__prefiltering_rename_vaf_to_af#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__prefiltering_rename_vaf_to_af#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__prefiltering_rename_vaf_to_af#


## prefiltering_artifact_homopolymer_annotation
Annotate variants with their frequency in the artifact panel (`Artifact`, `ArtifactMedian`, `ArtifactNrSD`).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__prefiltering__prefiltering_artifact_homopolymer_annotation#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__prefiltering_artifact_homopolymer_annotation#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__prefiltering_artifact_homopolymer_annotation#


# QC

Rules in `workflow/rules/qc_ont.smk`.

## [qc_ont_pycoqc](https://github.com/a-slide/pycoQC)
Report on the sequencing run from the MinKNOW sequencing summary.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__qc_ont__qc_ont_pycoqc#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__qc_ont__qc_ont_pycoqc#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__qc_ont_pycoqc#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__qc_ont_pycoqc#


## [qc_ont_sequali](https://github.com/rhpvorderman/sequali)
Quality report on the reads kept and removed by Filtlong.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__qc_ont__qc_ont_sequali#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__qc_ont__qc_ont_sequali#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__qc_ont_sequali#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__qc_ont_sequali#


## [qc_ont_mosdepth_overlap](https://github.com/brentp/mosdepth)
Mean coverage per target.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__qc_ont__qc_ont_mosdepth_overlap#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__qc_ont__qc_ont_mosdepth_overlap#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__qc_ont_mosdepth_overlap#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__qc_ont_mosdepth_overlap#


## qc_ont_mosdepth_merge
Merge the per-target coverage summaries into one CSV file.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__qc_ont__qc_ont_mosdepth_merge#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__qc_ont_mosdepth_merge#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__qc_ont_mosdepth_merge#


## qc_ont_yield_per_pool
Estimated read count and share of reads per amplicon in each pool.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__qc_ont__qc_ont_yield_per_pool#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__qc_ont__qc_ont_yield_per_pool#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__qc_ont_yield_per_pool#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__qc_ont_yield_per_pool#


## [qc_ont_picard_bed_to_interval_list](https://snakemake-wrappers.readthedocs.io/en/v5.0.1/wrappers/picard/bedtointervallist.html)
Convert `amplicons.bed` to a Picard interval list.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__qc_ont__qc_ont_picard_bed_to_interval_list#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__qc_ont_picard_bed_to_interval_list#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__qc_ont_picard_bed_to_interval_list#


## qc_ont_mosdepth_exons_multiqc
Write the coverage per exon as a MultiQC custom-content table.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__qc_ont__qc_ont_mosdepth_exons_multiqc#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__qc_ont__qc_ont_mosdepth_exons_multiqc#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__qc_ont_mosdepth_exons_multiqc#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__qc_ont_mosdepth_exons_multiqc#


## [qc_ont_mosdepth_overlap_timestep](https://github.com/brentp/mosdepth)
Coverage per target for each MinKNOW BAM batch (research).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__qc_ont__qc_ont_mosdepth_overlap_timestep#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__qc_ont__qc_ont_mosdepth_overlap_timestep#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__qc_ont_mosdepth_overlap_timestep#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__qc_ont_mosdepth_overlap_timestep#


## qc_ont_mosdepth_merge_timestep
Merge the coverage per batch (research).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__qc_ont__qc_ont_mosdepth_merge_timestep#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__qc_ont__qc_ont_mosdepth_merge_timestep#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__qc_ont_mosdepth_merge_timestep#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__qc_ont_mosdepth_merge_timestep#


## qc_ont_plot_yield_timestep
Plot the cumulative coverage per amplicon over time (research).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__qc_ont__qc_ont_plot_yield_timestep#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__qc_ont__qc_ont_plot_yield_timestep#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__qc_ont_plot_yield_timestep#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__qc_ont_plot_yield_timestep#


# Report

Rules in `workflow/rules/results_report.smk`.

## [results_report_xlsx](https://xlsxwriter.readthedocs.io/)
Write the Excel report with variants and coverage.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__results_report__results_report_xlsx#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__results_report__results_report_xlsx#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__results_report_xlsx#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__results_report_xlsx#


## results_report_copy_bed
Copy the BED files used in the analysis.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__results_report__results_report_copy_bed#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__results_report_copy_bed#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__results_report_copy_bed#


# Reference workflow

Rules in `workflow/rules/references_rules.smk`.

## [reference_rules_snv_indels_clairs_to_gvcf](https://github.com/HKU-BAL/ClairS-TO)
ClairS-TO with reference calls printed, for the background panel.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__references_rules__reference_rules_snv_indels_clairs_to_gvcf#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__references_rules__reference_rules_snv_indels_clairs_to_gvcf#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__reference_rules_snv_indels_clairs_to_gvcf#


## [reference_rules_snv_indels_deepsomatic_gvcf](https://github.com/google/deepsomatic)
DeepSomatic with gVCF output, for the background panel.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__references_rules__reference_rules_snv_indels_deepsomatic_gvcf#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__references_rules__reference_rules_snv_indels_deepsomatic_gvcf#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__reference_rules_snv_indels_deepsomatic_gvcf#


## [reference_rules_align_unfiltered](https://software-docs.nanoporetech.com/dorado/latest/basecaller/alignment/)
Align reads not filtered by Filtlong.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__references_rules__reference_rules_align_unfiltered#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__references_rules__reference_rules_align_unfiltered#


## [reference_rules_bam_sort_unfiltered](https://www.htslib.org/doc/samtools-sort.html)
Sort the unfiltered aligned reads.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__references_rules__reference_rules_bam_sort_unfiltered#


## [reference_rules_bam_index_unfiltered](https://www.htslib.org/doc/samtools-index.html)
Index the unfiltered sorted BAM file.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__references_rules__reference_rules_bam_index_unfiltered#


## create_background_file_longread_vardict_unfiltered
Background panel from VarDict gVCFs of unfiltered reads.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__references_rules__create_background_file_longread_vardict_unfiltered#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__references_rules__create_background_file_longread_vardict_unfiltered#


## create_artifact_file_longread
Artifact panel from the soft-filtered ensemble VCFs of the panel samples.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__references_rules__create_artifact_file_longread#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__references_rules__create_artifact_file_longread#

