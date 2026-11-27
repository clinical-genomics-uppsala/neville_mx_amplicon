# Excel report

`results_report_xlsx` (`workflow/scripts/results_report_xlsx.py`, tables built in `results_report_create_tables.py`) writes one
workbook per sample: `Results/Reports/{experiment}_{sample}_{type}.xlsx`.

## Input

| Input | Content |
|-------|---------|
| `*.ensembled.artifact_annotated.filter.somatic_soft.vcf.gz` | Soft-filtered SNVs and indels. |
| `*.sniffles2.bcftools_view.vep_annotated.vcf.gz` | Annotated structural variants. |
| Mosdepth outputs | Coverage per amplicon, per target and per exon. |
| `*_coverage_per_amplicon.csv`, `*_yield_pool_{pooln}.csv` | Estimated read counts per amplicon and per pool. |
| `config_soft_filter_somatic.yaml`, `config_hard_filter_somatic.yaml` | Filter descriptions printed in the SNV sheet. |

## Sheets

### Overview
Sample and run IDs, processing date, fields for sign-off (created by, signed by, document number), links to the other sheets,
the average coverage over `amplicons.bed`, and the BED files and artifact panel used.

### SNV_Indels
The filters applied, then one row per variant with AF above 0.5% and length up to 50 bp:

Filter, RunID, DNAnr, Gene, Chr, Pos, Ref, Alt, AF, DP, Transcript, Mutation cds, ENSP, Consequence, COSMIC ids on pos,
Clinical Significance, dbSNP, Max Pop AF, Max Pop, Artifact Medians (ClairS-TO; DeepSomatic; VarDict),
Artifact calls (ClairS-TO; DeepSomatic; VarDict), Callers.

The `Filter` column holds the soft-filter flags (see [Soft filtering](snv_indels.md#soft-filtering)); `PASS` variants have none.

### SVs
Variants from Sniffles2 with AF above 0.5%, and variants longer than 50 bp from the small-variant callers. Same columns as the
SNV sheet, with an extra `SV length` column and without the artifact columns.

### Coverage
* Average coverage per amplicon, including overlapping regions;
* Estimated read count per amplicon (non-overlapping regions);
* One table per pool with the read count and percentage of reads per amplicon (see [Yield per amplicon and per pool](qc.md#yield-per-amplicon-and-per-pool)).

### Exon Coverage
Average coverage per exon, from `exon_coverage.bed`.
