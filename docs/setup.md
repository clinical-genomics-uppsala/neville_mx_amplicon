# Pipeline setup and configuration

## Dependencies

| Software | Version |
|----------|---------|
| [Python](https://www.python.org/) | 3.9 or later |
| [Snakemake](https://snakemake.readthedocs.io/en/stable/) | 7.x (`>=7.8.0,<8`) |
| [hydra-genetics](https://hydra-genetics.readthedocs.io/en/latest/) | 3.0.0 |
| [Singularity](https://sylabs.io/docs/) / Apptainer | 3.x |
| [pandas](https://pandas.pydata.org/) | `>=1.3.1` |

The other Python packages are listed in `requirements.txt`.

!!! warning
    `smart-open` must be older than 7.0, otherwise the Picard rules fail.

All tools run in containers, so the only software needed on the host is Python, Snakemake and Singularity/Apptainer.
Basecalling with Dorado needs a GPU node with CUDA (see [GPU basecalling](running.md#gpu-basecalling)).

## Installation

```bash
git clone https://github.com/clinical-genomics-uppsala/neville_mx_amplicon.git
cd neville_mx_amplicon
python3 -m venv .venv
source .venv/bin/activate
pip install -I -r requirements.txt
mkdir tmp   # temporary files from some variant callers
```

## Configuration files

| File | Content |
|------|---------|
| `config/config.yaml` | Paths to references and BED files, amplicons and pools, run information, containers and tool parameters. |
| `config/resources.yaml` | Cluster resources (threads, memory, time, partition, GPU) per rule. |
| `config/output_files.yaml` | Files copied to `Results/` (see [Result files](result_files.md)). |
| `config/config_soft_filter_somatic.yaml` | Soft filters applied to the annotated ensemble VCF (see [SNVs and indels](snv_indels.md#soft-filtering)). |
| `config/config_hard_filter_somatic.yaml` | Hard filters, listed in the Excel report. |
| `config/rename_vaf_to_af.txt` | Rename map used to turn `FORMAT/VAF` into `FORMAT/AF`. |
| `config/multiqc_config.yaml` | MultiQC report layout. |
| `config/config_references.yaml`, `config/output_files_references.yaml`, `config/resources_references.yaml` | Settings for the [reference workflow](running_references.md). |

Both `config.yaml` and `resources.yaml` are validated against the schemas in `workflow/schemas/` when the pipeline starts.
The parameters available for each rule are listed in the [rule list](softwares.md).

### Paths to edit

| Key | Description |
|-----|-------------|
| `ref_data` | Reference genome FASTA (GRCh38) used by Dorado, ClairS-TO and DeepSomatic. |
| `bed_files` | Directory holding all the BED files of the design (see [Reference files](references.md#bed-files)). |
| `reference.fasta`, `reference.fai`, `reference.sequence_dict` | Reference genome, its index and its sequence dictionary. |
| `reference.design_bed` | BED file with the amplicons. |
| `reference.design_intervals` | Picard interval list of the amplicons. It is created from `amplicons.bed` if missing. |
| `reference.artifacts` | Artifact panel TSV built by the [reference workflow](running_references.md). |
| `vep.vep_cache` | Directory with the VEP cache. |
| `*.bed_file` | Caller-specific BED files (`snv_indels_svs_clairs_to`, `snv_indels_svs_deepsomatic`, `vardict`, `prefiltering_bcftools_include_region`). |

### Run information

These keys are usually given on the command line with `--config` (see [Running the pipeline](running.md)).

| Key | Description |
|-----|-------------|
| `runfolder` | Path to the folder holding the sequencing runs. |
| `batchid` | Batch/experiment ID; becomes the `{experiment}` wildcard. |
| `runid` | MinKNOW run folder name (e.g. `20251017_1440_MN48987_FBD39643_37b961fb`). |
| `multisample` | `True` for barcoded, multiplexed samples; `False` for one sample per flowcell. |
| `samplesheet` | Dorado sample sheet (CSV), needed when `multisample` is `True`. |
| `raw_data` | Name of the POD5 folder inside the run folder (default `pod5`). |
| `timestep_minknow` | Interval in minutes between the BAM batches written by MinKNOW, used for the yield-over-time plots. |

The POD5 input folder is `{runfolder}/{batchid}/{runid}/{raw_data}`.

### Amplicons and pools

```yaml
amplicons:
  ["FLT3_ITD_3kb_1", "FLT3_TKD_3kb_1", "IDH1_1", "IDH2_1", "NPM1_5",
   "TP53_3kb_A1", "TP53_3kb_B1", "TP53_C10", "TP53_J3", "TP53_D2", "TP53_E1", "TP53_F2", "TP53_G4", "TP53_H4", "TP53_I4"]

extra_regions:
  ["TP53_3kb_A1_only", "TP53_3kb_B1_only", "TP53_C10_only", "TP53_D2+J3", "TP53_D2_only",
   "TP53_E1_only", "TP53_F2_only", "TP53_G4_only", "TP53_H4_only", "TP53_I4_only"]

pools:
  1: ["TP53_3kb_A1_only", "TP53_D2_only", "TP53_F2_only", "TP53_I4_only", "FLT3_TKD_3kb_1"]
  2: ["IDH1_1", "TP53_3kb_B1_only", "TP53_G4_only", "TP53_D2+J3", "NPM1_5"]
  3: ["IDH2_1", "TP53_C10_only", "TP53_E1_only", "TP53_H4_only", "FLT3_ITD_3kb_1"]
```

* `amplicons`: one entry per amplicon. Each name must match a BED file `{bed_files}/{name}.bed`.
* `extra_regions`: sub-regions of overlapping *TP53* amplicons that are covered by one amplicon only (`*_only`),
  plus the region shared by D2 and J3 (`TP53_D2+J3`). Each one also needs its own BED file.
* `pools`: which amplicons (or sub-regions) belong to each PCR pool. Used to compute the yield per pool.

Change these three keys and the matching BED files to use another amplicon design.

### Containers

`default_container` is used by rules that do not set their own. Every other rule takes its image from `<rule_name>.container`.
The main images are:

| Tool | Version | Container |
|------|---------|-----------|
| Dorado | 1.1.1 | `docker://hydragenetics/dorado:1.1.1` |
| Samtools | 1.21 | `docker://hydragenetics/samtools:1.21` |
| Filtlong | 0.2.1 | `docker://hydragenetics/filtlong:0.2.1` |
| ClairS-TO | 0.4.0 | `docker://hkubal/clairs-to:v0.4.0` |
| DeepSomatic | 1.9.0 | `docker://google/deepsomatic:1.9.0` |
| VarDict | 1.8.3 | `docker://hydragenetics/vardict:1.8.3` |
| Sniffles2 | 2.6.1 | `docker://hydragenetics/sniffles2:2.6.1` |
| VT | 2015.11.10 | `docker://hydragenetics/vt:2015.11.10` |
| bcbio-variation-recall | 0.2.6 | `docker://hydragenetics/bcbio-vc:0.2.6` |
| VEP | 111.0 | `docker://hydragenetics/vep:111.0` |
| Mosdepth | 0.3.6 | `docker://hydragenetics/mosdepth:0.3.6` |
| Picard | 2.25.4 | `docker://hydragenetics/picard:2.25.4` |
| Sequali | 0.12.0 | `docker://hydragenetics/sequali:0.12.0` |
| pycoQC | 2.5.0.23 | `docker://hydragenetics/pycoqc:2.5.0.23` |
| MultiQC | 1.25.1 | `docker://hydragenetics/multiqc:1.25.1` |

!!! note
    DeepSomatic needs the bind mount `/usr/lib/locale/:/usr/lib/locale/` in the Singularity arguments.

## Sample files

The pipeline reads `samples.tsv` and `units.tsv` from the working directory (keys `samples_run` and `units_run`).
Create them with [`hydra-genetics create-input-files`](https://hydra-genetics.readthedocs.io/en/latest/run_pipeline/create_sample_files/):

```bash
hydra-genetics create-input-files -d <folder_with_merged_bam> -t T -p ONT -f
```

For ONT data, the command reads the run information from the header of a basecalled, **unaligned** BAM file.
Merge the MinKNOW BAM files first (`samtools merge`), and put the merged file in a folder with no sub-folders.

| Column | Description |
|--------|-------------|
| **`samples.tsv`** | |
| sample | Unique sample ID, one per row. |
| **`units.tsv`** | |
| sample | Same sample ID as in `samples.tsv`. |
| type | Sample type, `T` for this analysis. |
| platform | Sequencing platform, `ONT`. |
| machine | Machine ID. |
| processing_unit | Flowcell ID. |
| run_id | Sequencing run ID. |
| barcode | Library barcode, `NA` if none. |
| methylation | Whether methylation information is used, e.g. `No`. |
| basecalling_model | Basecalling model used **in MinKNOW**, e.g. `dna_r10.4.1_e8.2_400bps_sup@v5.0.0`. |
| bam | Path to the unaligned BAM file. |

Examples are in `.tests/integration/config/`.

## Profiles

Snakemake profiles hold the executor settings and the Singularity arguments.

| Profile | Use |
|---------|-----|
| `profiles/local/` | Local run, e.g. on a laptop. |
| `profiles/marvin_drmaa/` | CGU cluster (Slurm through DRMAA). |
| `profiles/miarka_slurm/` | Slurm cluster with GPU partition for basecalling. |

GPU rules need GPU resources from the scheduler and `--nv` in the Singularity arguments. With DRMAA, set `gres`
(e.g. `gres: "--gres=gres:gpu:4"`) and the GPU `partition` for the basecalling rules in `resources.yaml`; the DRMAA string in
`profiles/marvin_drmaa/config.yaml` passes `{resources.gres}` to Slurm. With the Slurm executor, use `set-resources` in the
profile as in `profiles/miarka_slurm/config.yaml`.
