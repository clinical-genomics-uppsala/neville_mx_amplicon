# Building reference panels

`workflow/Snakefile_references` builds panels from a set of reference samples. These panels describe the sequencing and
variant-calling noise of the assay:

* the **artifact panel** lists the variants called in the reference samples, with their frequency per caller.
  The main pipeline uses it to flag recurrent artifacts (see [Artifact annotation](snv_indels.md#artifact-annotation));
* the **background panels** give, for each position, the median and standard deviation of the non-reference allele
  frequency seen in the reference samples, one panel per caller. They are used for noise analyses; the main pipeline
  does not use them yet.

The reference workflow imports every rule of the main `Snakefile` (except `all`) as a module, so the same preprocessing,
alignment and variant calling are applied to the reference samples.

## Input

| Config key (`config/config_references.yaml`) | Description |
|------|-------------|
| `samples_ref`, `units_ref` | Sample files of the samples to process. |
| `panel_ref` | TSV with columns `experiment`, `sample`, `type`: the samples to include in the panels. |
| `output` | `config/output_files_references.yaml`; files are written to `reference_files/`. |
| `reference_rules_snv_indels_*_gvcf` | Containers, parameters and BED file (`amplicons-background.bed`) for the gVCF calling. |
| `create_background_file.min_dp` | Positions with depth at or below this value are ignored (default 300). |
| `create_background_file.max_af` | Allele frequencies above this value are treated as real variants and ignored (default 0.03). |

## Steps

| Rule | Description |
|------|-------------|
| `reference_rules_snv_indels_clairs_to_gvcf` | ClairS-TO with `--print_ref_calls`, to get allele counts at every position. |
| `reference_rules_snv_indels_deepsomatic_gvcf` | DeepSomatic with `--output_gvcf`. |
| `reference_rules_snv_indels_vardict_gvcf` | VarDict with `-p` (report all positions) on the filtered, soft-clipped BAM file. |
| `reference_rules_align_unfiltered`, `reference_rules_bam_sort_unfiltered`, `reference_rules_bam_index_unfiltered`, `reference_rules_snv_indels_vardict_gvcf_unfiltered` | Same VarDict gVCF on reads **not** filtered by Filtlong, to measure the effect of filtering on noise. |
| `references_create_background_file_{clairs_to,deepsomatic,vardict}` | Background panel per caller, from the hydra-genetics [references](https://hydra-genetics-references.readthedocs.io/en/latest/) module. |
| `create_background_file_longread_vardict_unfiltered` | Background panel from the unfiltered VarDict gVCFs (`workflow/scripts/create_background_file_longread.py`). |
| `create_artifact_file_longread` | Artifact panel from the soft-filtered ensemble VCFs (`workflow/scripts/create_artifact_file_longread.py`). |

### Artifact panel

For each variant (chromosome, position, type, ref, alt) in the ensemble VCFs of the panel samples, and for each caller listed
in the `CALLERS` INFO field, the script counts the observations and collects the allele frequencies (`FORMAT/AF`).
The median and standard deviation are computed only when a variant has at least 4 observations; otherwise they are set to 0 and 1000.

### Background panel

For each position in the gVCFs with a depth above `min_dp`, the non-reference allele frequency is computed from `FORMAT/AD`.
Frequencies above `max_af` (or above `1 - max_af`, for homozygous positions) are skipped as likely real variants.
The panel reports the median and standard deviation per position, or the mean with SD 0 when there are fewer than 4 observations.

## Running

`workflow/scripts/start_marvin_multisample_references.sh` prepares the sample files as for the main pipeline, then runs:

```bash
snakemake --profile profiles/marvin_drmaa/ -s workflow/Snakefile_references \
    --configfiles config/config.yaml config/config_references.yaml \
    --config runfolder=<runFolder> batchid=<batchId> runid=<runId> multisample=True samplesheet=<samplesheet.csv> \
    --notemp --latency-wait 30
```

!!! note
    The background panel rules read the per-sample gVCFs from `reference_files/`, which is where the first pass of the
    workflow copies them. Panels covering several runs are therefore built in two steps: run the workflow on each run to
    collect the gVCFs, then list all the samples in `panel_ref` and run it again.

## Output (`reference_files/`)

| File | Description |
|------|-------------|
| `artifact_panel.tsv` | Artifact panel; set `reference.artifacts` in `config/config.yaml` to this file. |
| `background_panel_{clairs_to,deepsomatic,vardict}.tsv` | Background panel per caller. |
| `background_panel_vardict_unfiltered.tsv` | Background panel from unfiltered reads. |
| `{experiment}_{sample}_{type}_snv.clairs_to.g.vcf.gz`, `.deepsomatic.g.vcf.gz`, `.vardict.g.vcf.gz`, `.unfiltered_vardict.g.vcf.gz` | gVCFs per sample. |
| `{experiment}_{sample}_{type}.ont_reads.snv-indels.vcf.gz` | Ensemble VCF per sample. |

## Noise analysis scripts

Two command-line scripts analyse the panels and the BAM files at known hotspot positions. They are not part of the Snakemake workflow.
Arguments can be given in a file with `@args.txt`; templates are in `workflow/scripts/*_args.template.txt`.

| Script | Description |
|--------|-------------|
| `references_report_noise_artifacts.py` | For each hotspot in a CSV file, collects the background and artifact panel values per caller and writes them to CSV and JSON. |
| `hotspots_seq_noise_analysis.py` | For each hotspot, counts the reads carrying the reference and alternative allele in each BAM file of the reference panel and computes the VAF. |
