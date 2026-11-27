# Result files

The pipeline generates many intermediate files, most of them temporary. The files kept after the run are listed in
`config/output_files.yaml` and copied to `Results/` by automatically generated `_copy_*` rules.
Add or remove entries in that file to change what is delivered.

`{experiment}` is the batch ID of the run, `{sample}` the sample ID and `{type}` the sample type (`T`).

## Data

| File | Description |
|------|-------------|
| `Results/Data/{experiment}_{sample}_{type}.ont_reads.unfiltered.fastq.gz` | Basecalled, adapter-trimmed reads before filtering. |
| `Results/Data/{experiment}_{sample}_{type}.ont_reads.filtered.fastq.gz` | Reads kept by Filtlong (quality and length filtering). |
| `Results/Data/{experiment}_{sample}_{type}.ont_reads.filtered.aligned.sorted.soft-clipped.bam` | Aligned, sorted reads with primers soft-clipped. |
| `Results/Data/{experiment}_{sample}_{type}.ont_reads.filtered.aligned.sorted.soft-clipped.bam.bai` | BAM index. |
| `Results/Data/{experiment}_{sample}_{type}.ont_reads.snv-indels.vcf` | Ensemble SNVs and indels from ClairS-TO, DeepSomatic and VarDict, VEP- and artifact-annotated, soft-filtered. |
| `Results/Data/{experiment}_{sample}_{type}.ont_reads.cnv-sv.vcf.gz` | Structural variants from Sniffles2, restricted to the target regions. |

## Reports

| File | Description |
|------|-------------|
| `Results/Reports/{experiment}_{sample}_{type}.xlsx` | Excel report with variants and coverage tables. See [Excel report](report.md). |
| `Results/Reports/{experiment}_{sample}_{type}_multiqc_amplicons.html` | MultiQC report with read quality, coverage, on-target and alignment metrics. |
| `Results/Reports/{experiment}_{sample}_{type}.ont_reads.filtered.fastq.gz.html` | Sequali read quality report for the filtered reads. |
| `Results/Reports/{experiment}_{sample}_{type}_sequencing_report_pycoqc.html` | pycoQC report of the sequencing run, from the MinKNOW `sequencing_summary*.txt`. |

## Logs and benchmarks

Log files (`*.log`) and benchmark files (`*.benchmark.tsv`) are written next to each rule's output in the working directory
(`basecalling/`, `prealignment/`, `alignment/`, `snv_indels/`, `cnv_sv/`, `results/`, `reports/`) and are not removed.

## Reference workflow

The [reference workflow](running_references.md) uses `config/output_files_references.yaml` and writes its files to `reference_files/`.
See that page for the list of files.
