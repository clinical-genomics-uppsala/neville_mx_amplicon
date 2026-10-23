# Reference files

## Reference genome

| Config key | Description |
|------------|-------------|
| `ref_data`, `reference.fasta` | GRCh38 FASTA, e.g. `GCA_000001405.15_GRCh38_no_alt_analysis_set` (CGU uses a masked version). |
| `reference.fai` | FASTA index (`samtools faidx`). |
| `reference.sequence_dict` | Sequence dictionary (`picard CreateSequenceDictionary`), needed to build the Picard interval list. |

## BED files

All BED files live in the folder given by `bed_files`. Examples for the current design are in
`.tests/integration/test_data/bedfiles/`.

| File | Used by | Description |
|------|---------|-------------|
| `amplicons.bed` | Mosdepth, Picard, Excel report | One line per amplicon. **Must be sorted by chromosomal order, with no empty line at the end.** |
| `amplicons-primers.bed` | `samtools ampliconclip` | Position of the forward and reverse primer of each amplicon; these bases are soft-clipped. |
| `amplicons.interval_list` | Picard HsMetrics | Picard interval list; created from `amplicons.bed` if it does not exist (`reference.design_intervals`). |
| `{target}.bed` | Split BAM, Mosdepth | One file per entry in `amplicons` and `extra_regions`. For overlapping *TP53* amplicons, `*_only.bed` holds the part covered by that amplicon only, and `TP53_D2+J3.bed` the region shared by D2 and J3. |
| `amplicons-clairs-to.bed` | ClairS-TO | Regions where ClairS-TO calls variants. |
| `amplicons-exons-deepsomatic.bed` | DeepSomatic | Regions where DeepSomatic calls variants (exons of the targets). |
| `NPM1-FLT3-vardict.bed` / `amplicons-vardict.bed` | VarDict | VarDict regions in amplicon mode. Needs extra columns, see `vardict.bed_columns`. |
| `amplicons-sniffles2.bed` | bcftools | Sniffles2 cannot be restricted to regions, so its calls are filtered on these regions afterwards. |
| `exon_coverage.bed` | Mosdepth | Exons for which the mean coverage is reported (`mosdepth_exons.bed`). |
| `amplicons-background.bed` | Reference workflow | Regions where gVCFs are computed for the background panel. |

## Dorado models

Basecalling models (e.g. `dna_r10.4.1_e8.2_400bps_sup@v5.2.0`) are downloaded by `dorado download` if missing. For multiplexed runs,
the model is first looked up in the folder given by `DORADO_MODELS`, which is useful on nodes without internet access.
`dorado_model_extra` names an extra model folder copied together with the basecalling model.

## VEP cache

`vep.vep_cache` points at the VEP cache (`--offline --cache --merged` by default). The cache must match the VEP version (111) and GRCh38.

## Artifact panel

`reference.artifacts` is a TSV file with, for each position and variant seen in a set of reference samples, the number of
observations and the median and standard deviation of the allele frequency, per caller:

```
Chromosome  pos  variant_type  ref  alt  clairs_to  median_MAF  sd_MAF  deepsomatic  median_MAF  sd_MAF  vardict  median_MAF  sd_MAF
```

It is built with the [reference workflow](running_references.md) (`create_artifact_file_longread`) and used to annotate the
ensemble VCF (see [SNVs and indels](snv_indels.md#artifact-annotation)).
