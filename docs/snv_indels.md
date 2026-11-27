# SNVs and indels

Small variants are called with three callers, merged into one ensemble VCF, then annotated and soft-filtered.
See the hydra-genetics [snv_indels](https://hydra-genetics-snv-indels.readthedocs.io/en/latest/),
[annotation](https://hydra-genetics-annotation.readthedocs.io/en/latest/) and
[filtering](https://hydra-genetics-filtering.readthedocs.io/en/latest/) modules for details on the imported rules.

## Pipeline output files

* `Results/Data/{experiment}_{sample}_{type}.ont_reads.snv-indels.vcf`
* SNV_Indels sheet of the [Excel report](report.md)

## Variant calling

All callers run on the soft-clipped BAM file, in tumour-only mode, restricted to their own BED file.

### ClairS-TO
**[ClairS-TO](https://github.com/HKU-BAL/ClairS-TO)** v0.4.0 (`snv_indels_svs_clairs_to`) is a deep-learning somatic caller for long reads.
It writes separate SNV and indel VCF files, which `snv_indels_svs_clairs_to_concat` concatenates and sorts with bcftools.
The verdict module (germline/somatic/subclonal classification) is disabled with `--disable_verdict`.

| Option | Value | Description |
|--------|-------|-------------|
| `platform` | `ont_r10_dorado_sup_5khz` | Model matching R10.4.1 flowcells and Dorado SUP basecalling. |
| `snv_min_af` | 0.005 | Minimum AF for SNV candidates. |
| `indel_min_af` | 0.005 | Minimum AF for indel candidates. |
| `bed_file` | `amplicons-clairs-to.bed` | Regions to call. |

### DeepSomatic
**[DeepSomatic](https://github.com/google/deepsomatic)** v1.9.0 (`snv_indels_svs_deepsomatic`) runs with the `ONT_TUMOR_ONLY` model.

| Option | Value | Description |
|--------|-------|-------------|
| `extra` | `--make_examples_extra_args="vsc_min_fraction_snps=0.01,vsc_min_fraction_indels=0.01,small_model_snp_gq_threshold=20,small_model_indel_gq_threshold=20"` | Lower candidate AF thresholds to 1%. |
| `bed_file` | `amplicons-exons-deepsomatic.bed` | Regions to call (exons). |

DeepSomatic downsamples the reads before calling, so `AD` and `DP` in its VCF refer to the downsampled reads.

### VarDict
**[VarDict](https://github.com/AstraZeneca-NGS/VarDict)** v1.8.3 (`snv_indels_vardict`, from hydra-genetics) was designed for
short reads, but its amplicon mode (`-a`) works with ONT data. It gives many false-positive short indels at low VAF.

| Option | Value | Description |
|--------|-------|-------------|
| `allele_frequency_threshold` | 0.002 | Minimum AF reported. |
| `extra` | `-k 0 -L 1000 -a 50:0.95` | No local realignment; amplicon mode. |
| `bed_file` | `NPM1-FLT3-vardict.bed` | Regions to call, with the extra columns needed in amplicon mode. |
| `bed_columns` | `-c 1 -S 2 -E 3 -g 4` | Column layout of the BED file. |

!!! note
    GATK Mutect2 (`snv_indels_gatk_mutect2`) is defined but not part of the ensemble. SAVANA (`snv_indels_svs_savana`) is kept
    for later use; it does not handle amplicon data yet ([savana#97](https://github.com/cortes-ciriano-lab/savana/issues/97)).

## Pre-filtering

`prefiltering_bcftools_view` removes the calls that each caller marks as false positives, before merging:

| Caller | bcftools expression |
|--------|---------------------|
| DeepSomatic | `-i 'QUAL>0'` |
| ClairS-TO | `-e 'FILTER=="LowQual" & QUAL==0'` |
| VarDict | `-i 'QUAL>120'` |

Set the expressions under `prefiltering_bcftools_view.filter.<caller>.extra`.

## Decomposition and normalisation

Multi-allelic variants are split with **[vt decompose](https://genome.sph.umich.edu/wiki/Vt#Decompose)** and indels are
left-aligned with **[vt normalize](https://genome.sph.umich.edu/wiki/Vt#Normalization)**, so that the callers' outputs can be compared.

## Ensemble calling

**[bcbio-variation-recall ensemble](https://github.com/bcbio/bcbio.variation.recall)** merges the three VCF files
(`snv_indels_bcbio_variation_recall_ensemble`).

| Option | Value | Description |
|--------|-------|-------------|
| `callers` | `clairs_to`, `deepsomatic`, `vardict` | Callers to merge, in order of priority. The soft filter on artifacts assumes this order (indices 0, 1, 2). |
| `support` | 1 | Number of callers needed to keep a variant. |

The `CALLERS` INFO field lists the callers that found each variant.

## Annotation

The ensemble VCF is annotated with **[VEP](https://www.ensembl.org/info/docs/tools/vep/index.html)** v111 (`annotation_vep`):

```
--offline --cache --merged --format vcf --assembly GRCh38 --check_existing --pick --variant_class --everything
--pick_order mane_select,mane_plus_clinical,canonical,biotype,rank,appris,tsl,ccds,length,ensembl,refseq
```

`prefiltering_rename_vaf_to_af` then renames `FORMAT/VAF` to `FORMAT/AF` with `bcftools annotate --rename-annots`
(rename map in `config/rename_vaf_to_af.txt`), so all records use `AF`.

## Artifact annotation

`prefiltering_artifact_homopolymer_annotation` (`workflow/scripts/artifact_and_homopolymer_annotation.py`) looks up each variant in
the [artifact panel](references.md#artifact-panel) and adds three INFO fields, with one value per caller in the order of the panel:

| Field | Description |
|-------|-------------|
| `Artifact` | Number of panel samples in which the variant was called, per caller, followed by the panel size. |
| `ArtifactMedian` | Median AF of the variant in the panel samples. |
| `ArtifactNrSD` | Number of standard deviations between the sample AF and the panel median: `|AF - median| / sd`. |

SNVs are matched on the exact substitution (e.g. `C>T`) and then on the generic type `SNV`; indels on `insertion`/`deletion`, then `INDEL`.
Variants absent from the panel get 0 observations and `ArtifactNrSD=1000`.

## Soft filtering

`filtering_filter_vcf` (hydra-genetics [filtering](https://hydra-genetics-filtering.readthedocs.io/en/latest/)) applies the filters in
`config/config_soft_filter_somatic.yaml`. Variants are kept, and the flag is added to the `FILTER` column:

| Flag | Condition |
|------|-----------|
| `Intron` | Intronic, unless also splice or with a COSMIC ID at the position. |
| `VAF` | `FORMAT/AF` below 0.01. |
| `Consequence` | Intergenic, NMD transcript, non-coding transcript, upstream/downstream, TFBS or regulatory region variant. |
| `PopAF_0.02` | Maximum population AF (1000 Genomes, ESP, gnomAD exomes) above 2%. |
| `Biotype` | Not annotated as `protein_coding`. |
| `Artifact_gt_5` | For any caller: seen in more than 5 panel samples and less than 5 SD from the panel median. |

The hard filters in `config/config_hard_filter_somatic.yaml` (gnomAD exome AF > 0.5%, AD < 50, not an SNV) are listed in
the Excel report.

## Allele frequency per caller

The callers compute the allele frequency differently.

**ClairS-TO**: `AF = AD_alt / DP`. For example, `GT:GQ:DP:AF:AD ... 0/1:25:5137:0.1468:4371,754` gives 754 / 5137 = 0.1468.

**DeepSomatic**: `AF = AD[1] / DP`, computed on the downsampled reads
([details](https://github.com/google/deepvariant/issues/880#issuecomment-2338514537)).

**Sniffles2** (structural variants): `VAF = SUPPORT / COVERAGE`, see [Structural variants](cnv_sv.md#allele-frequency).
