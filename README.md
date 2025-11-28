# <img src="images/hydragenetics.png" width=40 /> neville_mx_amplicon

#### A pipeline to filter, align, and analyze Nanopore sequence data from pooled amplicons

![Lint](https://github.com/hydra-genetics/pipeline_pool_amplicon/actions/workflows/lint.yaml/badge.svg?branch=develop)
![Snakefmt](https://github.com/hydra-genetics/pipeline_pool_amplicon/actions/workflows/snakefmt.yaml/badge.svg?branch=develop)
![snakemake dry run](https://github.com/hydra-genetics/pipeline_pool_amplicon/actions/workflows/snakemake-dry-run.yaml/badge.svg?branch=develop)
![integration test](https://github.com/hydra-genetics/pipeline_pool_amplicon/actions/workflows/integration.yaml/badge.svg?branch=develop)

![pycodestyle](https://github.com/hydra-genetics/pipeline_pool_amplicon/actions/workflows/pycodestyle.yaml/badge.svg?branch=develop)
![pytest](https://github.com/hydra-genetics/pipeline_pool_amplicon/actions/workflows/pytest.yaml/badge.svg?branch=develop)

[![License: GPL-3](https://img.shields.io/badge/License-GPL3-yellow.svg)](https://opensource.org/licenses/gpl-3.0.html)

/!\ README under construction

## :speech_balloon: Introduction

This module implements a workflow with Snakemake for the analysis of Nanopore sequence data from pooled amplicons
at Clinical Genomics Uppsala (CGU) in view of supporting the diagnostics of Acute Myeloid Leukemia (AML).
The module consists of ... build upon [hydra-genetics]().
Multiplexed amplicon sequencing: severral targets are amplified in **multiplexed PCR** settings. 
The **sample** to be analyzed is **not multiplexed** with other ones.

Aim with the multiplexed setup: decrease the workload of the amplification step in the wet lab.
Caution: if two amplicons cover regions that overlap, those amplicons must be assigned to different pools.

The data were obtained from sequencing the prepared library on a MinION machine with the sequencing kit 14.
The duplex basecalling is performed with [Dorado](https://github.com/nanoporetech/dorado) 
([documentation](https://dorado-docs.readthedocs.io/en/latest/)).
Dorado is the currently recommended software for MinION output.

### Steps performed in the analysis

The main processing steps are:

1. Basecalling with *Dorado* from POD5 files (raw output of ONT sequencing machines),
> You may skip this step and directly use the basecalled BAM files provided by MinKNOW in the `bam_pass` directory
> of the sequencing run folder. 
> In this case, you need to edit the rules in `workflow/rules/basecalling.smk` accordingly.
> You may first use `samtools merge <output.bam> <input1.bam> <input2.> bam> <...>` 
> to merge the BAM files in the `bam_pass` directory into a single BAM file.

2. As only one sample is sequenced at a time, the data do not carry any barcode neither they need to be demultiplexed.
The "demultiplexing" of the amplicons is done later after aligning the reads, based on the coordinates of the targets.

3. Filtering with *Filtlong*: 
only reads with sufficient quality (based on the Phred quality scores) are kept. No reference is used. 
Reads that are too short or too long are filtered out.
Very long reads (over 4kb in our case) can be chimeras that result from concatenated amplicons.
The `--split` option could be used as alternative strategy, in such case a reference genome is required.
The thresholds for minimal and maximal length may be adjusted depending on the set of amplicons.
See the [documentation for Filtlong](https://github.com/rrwick/Filtlong/tree/main) for more options.

4. Alignment with *Dorado* and soft-clipping with *samtools*,

5.Variant calling with:
  - the software [ClairS-TO](https://github.com/HKU-BAL/ClairS-TO) via the Docker container provided by the development team,
  - the software [DeepSomatic](https://github.com/google/deepsomatic) via Docker container provided by the development team and the [documented examples](https://github.com/google/deepsomatic/blob/r1.8/docs/deepsomatic-case-study-ont-tumor-only.md),
  - the software VarDict [](),
  - the software Sniffles v2 [](),

6. Decomposition of the variants with the software _VT_ and annotation of the variants with _VEP_,

7. Computation of some quality metrics and report them in an Excel file as well as in a HTML file with MultiQC:
  - Read quality with Sequali,
  - Estimated number of reads for each amplicon: approximated by the mean coverage computed with *mosdepth*,
  - Counts and proportion of reads for each amplicon relatively to the pool it belongs to,

8. Reporting the results of variant calling and variant filtering in an Excel file.

Estimating the number of reads for each amplicon is desirable for two reasons.
Firstly, we want to assess within each pool whether the amplified reads are well-balanced between the amplicons. 
For example, if a pool involves 5 targets, each amplicon is expected to represent ca. 20% of the reads in this pool.
Whenever the percentage of reads for an amplicon is much higher than the expected share, it indicates that this amplicon,
for some reason, takes over during the multiplexed PCR.
Conversely, an amplicon might be disfavored in the PCR stage against the other targets and results in very few amplified reads.
Secondly, we want to control the balance between the pools. That is, each of the 3 pools should yield after multiplexed PCR 
about 1/3 ~ 33% of the total number of reads. If not, one may want to readjust the quantities pipetted for each pool previously to sequencing.

Note that estimating the read counts from mean coverage is tricky in regions with overlapping amplicons as TP53.
Subregions covered by only one amplicon should be identified previously to execute mosdepth analysis in those subregions.

### Remark about the variant allele frequency calculated by ClairS-TO

In the [documentation for ClairS](https://github.com/HKU-BAL/ClairS/blob/main/docs/verdict.md):
* "AFsomatic is the expected allele frequency of the variant being somatic, calculated as p * V / (p * C + 2 * (1-p)), where p is the tumor purity, C is the copy number, and V is the variant allele count in the tumor"
* "AFgermline is the expected allele frequency of the variant being a germline, calculated as (p * V + (1-p))/ (p * C + 2 * (1-p))"

Examples of calculation of AF with variants tagged as somatic by ClairS-TO: 
1. `chr17   7670685 .       G       A       25.5208 PASS    FAU=294;FCU=0;FGU=1809;FTU=0;RAU=460;RCU=0;RGU=2562;RTU=2;SB=0.32822    GT:GQ:DP:AF:AD:AU:CU:GU:TU      0/1:25:5137:0.1468:4371,754:754:0:4371:2`
AF = AD_alt / DP = 754 / 5137 = 0.1468

2. `chr17   7675199 .       G       A       33.6370 PASS    FAU=3315;FCU=0;FGU=1048;FTU=0;RAU=1145;RCU=1;RGU=395;RTU=0;SB=0.30566   GT:GQ:DP:AF:AD:AU:CU:GU:TU      0/1:33:5910:0.7547:1443,4460:4460:1:1443:0`
AF = AD_alt / DP = 4460 / 5910 = 0.7547

## :heavy_exclamation_mark: Dependencies

In order to use this module, the following dependencies are required:

[![hydra-genetics](https://img.shields.io/badge/hydragenetics-v1.3.0-blue)](https://github.com/hydra-genetics/)
[![pandas](https://img.shields.io/badge/pandas-2.2.2-blue)](https://pandas.pydata.org/)
[![python](https://img.shields.io/badge/python-3.9-blue)]()
[![snakemake](https://img.shields.io/badge/snakemake-7.32.4-blue)](https://snakemake.readthedocs.io/en/stable/)
[![singularity](https://img.shields.io/badge/singularity-3.0.0-blue)](https://sylabs.io/docs/)

Download dorado_models locally.

The Python package `smart-open` must be version < 7.0, otherwise the analyses with Picard fail.
Idem with the package `snakemake-wrapper-utils`

Motivation for using Picard tools on long-read data: https://www.agilent.com/cs/library/applications/an-long-read-sureselect-xt-hs2-5994-7612en-agilent.pdf

## :school_satchel: Preparations

### Set of amplicons and pools

| Pool | Amplicon name | Amplicon length | Chromosome |
| --- | --- | --- | --- |
| 1 | TP53_D2 | 3000 | chr17 |

### BED files

A BED file with the physical position of each amplicon, **the amplicons must be sorted by chromosomal order and no empty line at the end of the file**. 
A BED file with the physical position of the forward and the reverse primer for each amplicon 
(required for soft-clipping with *samtools*).

In addition, or each overlapping amplicon, a BED file with the coordinates of the region 
covered by this amplicon **only** is necessary.

Caller-specific BED files: ddepsomatic, vardict, sniffles to

See examples of BED files in ...

### Sample data

Files created with hydra-genetics from the BAM-files generated the MinKNOW software every 10th minute:
```bash
samtools merge ...
```

Input data should be added to [`samples.tsv`](https://github.com/hydra-genetics/pipeline_pool_amplicon/blob/develop/config/samples.tsv)
and [`units.tsv`](https://github.com/hydra-genetics/pipeline_pool_amplicon/blob/develop/config/units.tsv).
The following information need to be added to these files:

| Column Id | Description |
| --- | --- |
| **`samples.tsv`** |
| sample | unique sample/patient id, one per row |
| **`units.tsv`** |
| sample | same sample/patient id as in `samples.tsv` |
| type | data type identifier (one letter), can be one of **T**umor, **N**ormal, **R**NA |
| platform | type of sequencing platform, e.g. `NovaSeq` |
| machine | specific machine id, e.g. NovaSeq instruments have `@Axxxxx` |
| flowcell | identifer of flowcell used |
| lane | flowcell lane number |
| barcode | sequence library barcode/index, connect forward and reverse indices by `+`, e.g. `ATGC+ATGC` |
| fastq1/2 | absolute path to forward and reverse reads |
| adapter | adapter sequences to be trimmed, separated by comma |

### Configuration of the pipeline

Edit the file paths and names in `config/config.yaml`.

Create a `tmp` directory in the root of the repository to store temporary files:
```bash
cd <root_of_the_repository>
mkdir tmp
```

### Set up the virtual environment

Create a virtual environment and install the required Python packages:
```bash
cd <root_of_the_repository>
python3 -m venv .venv
source .venv/bin/activate
pip install -I -r requirements.txt
```

### Set up the profile to run the pipeline

Modify the file `profiles/slurm/config.yaml` to match your configuration for the run.
At CGU, we use a cluster with Slurm installed on it and Slurm is interfaced with the application DRMAA.

The file `profiles/local/config.yaml` shows an example of profile that is used for local settings, for instance on a 
laptop.

**Note:** If you want to run the pipeline with basecalling executed on a GPU node, 
you must provide the correct resources to use the GPU e.g. `gres: "--gres=gres:gpu:4"` and pass them to DRMAA.
Moreover, you must pass the option `--nv` as an argument to Singularity/Apptainer. If not, you will encounter the 
following error:
error:
```
[2025-08-04 09:02:49.732] [info] Failed to load NVML
[2025-08-04 09:02:49.732] [error] device string set to cuda:all but no CUDA devices available.
CUDA device string format: "cuda:0,...,N" or "cuda:all".
terminate called after throwing an instance of 'std::runtime_error'
  what():  Could not open file: -
```
You may get the following message about no nv host found upon execution of the rules and allocation the resources 
for rules in which no GPU is involved:
```
INFO:    Could not find any nv files on this host!
```
which is not a problem as long as the GPU resources are correctly set up in the profile.

## :white_check_mark: Testing

### Linting
Run the linter to check the code style and the correctness of the Snakefile:
```bash
cd .tests/integration
snakemake --lint -n -s ../../workflow/Snakefile --configfiles config/config.yaml
```

### Dry run
Configure the virtual environment `.venv` and test the pipeline with a dry run to see if the rules are correctly 
defined and the dependencies are satisfied:

```bash
source .venv/bin/activate
cd .tests/integration
snakemake -n -s ../../workflow/Snakefile --configfile config/config.yaml --config runfolder=../../.tests/integration/test_data
```

**NB**: After creating a suitable [profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles), 
add the option `--profile profiles/slurm/` to run the pipeline on a cluster with SLURM installed on it.

### Small integration test
The small integration test may also serve as example of how to run the pipeline 
and how the working directory should be structured.

#### Test data and structure
The test sample D25-test007 is derived from a sample with known variants in the genes IDH1, IDH2, NPM1, FLT3, and TP53, 
which is a [gDNA mix](https://horizondiscovery.com/en/reference-standards/products/myeloid-gdna-reference-standard) of cell lines and synthetic DNA fragments that are commercially available from Horizon Discovery.
Disclaimer: the gDNA mix is not optimized for ONT sequencing and the results of the analysis are not sensible in a clinical way.
The point of the small integration test is only to demonstrate that the pipeline bioinformatically works.

The repository contains a setup for a small test in `.tests/integration` (structure shown below before running the test):

```bash
.tests/integration/
├── basecalling
│   └── dorado_duplex
│       └── D25-test007_T_reads.ont_adapt_trim.bam
├── config
│   ├── config_filter_bcftools.yaml
│   ├── config_hard_filter_somatic.yaml
│   ├── config_soft_filter_somatic.yaml
│   ├── config.yaml
│   ├── multiqc_config.yaml
│   ├── output_files.yaml
│   ├── rename_vaf_to_af.txt
│   ├── resources.yaml
│   ├── samples.tsv
│   └── units.tsv
├── dorado_models
│   └── dna_r10.4.1_e8.2_400bps_sup@v5.0.0
├── reference
│   ├── D25-test007_T_reads.ont_adapt_trim.filtered.aligned.bam
│   ├── TP53_chr17_GRCh38.dict
│   ├── TP53_chr17_GRCh38.fasta.gz
│   ├── TP53_chr17_GRCh38.fasta.fai
│   └── vep_cache
├── samples.tsv
├── test_data
│   ├── bedfiles
│   ├── D25-test007
│       └── ABC123
│           ├── bam_pass
│               ├── ABC123_pass_456_789_0.bam
│               ├── ABC123_pass_456_789_0.bam.bai
│               ├── ABC123_pass_456_789_1.bam
│               └── ABC123_pass_456_789_1.bam.bai
│           └── pod5
│   └── preprocessed
│       ├── D25-test007_T_reads.ont_adapt_trim.filtered.aligned.bam
│       └── D25-test007_T.ensembled.vep_annotated.vcf.gz
├── tmp
│   └── .
└── units.tsv
```

The folder `tmp` in `.tests/integration` must be created beforehand manually to store temporary files from some variant callers.

The test data are located in `.tests/integration/test_data/` and consist of:
* files structured as in a real sequencing run folder generated by MinKNOW, containing:
  * the BED files with the amplicon coordinates and primer coordinates,
  * the BAM files generated by MinKNOW in the `bam_pass` directory of the sequencing run folder.
  The file name *"ABC123_pass_456_789"* was chosen to match the usual ONT's naming for the basecalled files generated by MinKNOW.
  * the `pod5` directory is left empty as the POD5 files are too large to be stored on GitHub,
* preprocessed files for skipping the steps that do not need to be tested:
  * the basecalled BAM file located in `.tests/integration/basecalling/dorado_duplex/D25-test007_T_reads.ont_adapt_trim.bam`,
  with ca. 15,000 reads from amplicon sequencing of IDH1, IDH2, NPM1, FLT3, TP53 targets on a MinION device.
  This BAM file is used as input to the test pipeline, thus skipping the basecalling step.
  * the aligned BAM file located in `.tests/integration/test_data/preprocessed/D25-test007_T_reads.ont_adapt_trim.filtered.aligned.bam`
  that is used to overwrite the aligned BAM file generated by the alignment step in the test pipeline.
  This aligned BAM file contains exactly the same reads as in the input basecalled BAM file,
  but the reads were preliminarily aligned against the whole reference genome hg38.
  * VEP-annotated VCF files `.tests/integration/test_data/preprocessed/D25-test007_T_reads.*.vcf.gz`
  after variants were called with ClairS-TO, DeepSomatic, VarDict, and Sniffles2.
  This file is copied as output instead of running the VEP-annotation step in the test pipeline.
  
#### Preparation
1. Make sure that the virtual environment `.venv` is correctly set up as described above.

2. Verify that the input files `samples.tsv` and `units.tsv` are present in `.tests/integration/`, 
otherwise copy there the tsv files provided in `config`: `cp .tests/integration/config/*.tsv .tests/integration/`.

3. Set up the preprocessed basecalled BAM file in `.tests/integration/basecalling/dorado_duplex/`:
```bash
mkdir -p .tests/integration/basecalling/dorado_duplex/ && \
cp .tests/integration/test_data/preprocessed/D25-test007_T_reads.ont_adapt_trim.bam .tests/integration/basecalling/dorado_duplex/
```

4. If missing, create the temporary directory `.tests/integration/tmp/` to store temporary files: 
```bash
mkdir -p .tests/integration/tmp/
```

5. Check the connection to the internet. The test must be run in an environment which can connect to the internet. 
This is necessary to pull the Docker images into Singularity containers.

6. Decompress the reference genome if not done yet:
```bash
bgzip -d .tests/integration/reference/TP53_chr17_GRCh38.fasta.gz
```

#### Execution
The small integration test can be run as follows on a Linux-based OS
(modify the commands and the paths accordingly to your OS):

```bash
$ cd .tests/integration/
$ source ../../.venv/bin/activate
$ snakemake -s ../../workflow/Snakefile -j 1 --show-failed-logs --configfiles ../../config/config.yaml config/config.yaml  --use-singularity --singularity-args  " --cleanenv --containall --bind $PWD/tmp:/tmp -B $PWD:$PWD  -B $HOME -B /usr/lib/locale/:/usr/lib/locale/ --disable-cache "
```

The `--singularity-args` may be replaced by the ones that are suitable for your local OS if needed.
Optionally, you may add the option `--notemp` to keep all temporary files 
that are created during the run and inspect them afterwards. 
These options can be provided in the profile configuration file as well, this profile is OS-specific.

The whole execution takes ca. 2 minutes on a laptop with the following characteristics:
 * 22 CPU cores and 32 GiB RAM,
 * Processor: Intel® Core™ Ultra 7 165H × 22,
 * OS (64-bit): Ubuntu 22.04.5 LTS.

**Expected DAG output of the execution command**

```bash
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 1 (use --cores to define parallelism)
Rules claiming more threads will be scaled down.
Conda environments: ignored
Job stats:
job                                                             count
------------------------------------------------------------  -------
_copy_bed_files_amplicon                                            1
_copy_config_yaml                                                   1
_copy_excel_report                                                  1
_copy_filtered_fastq                                                1
_copy_mosdepth_coverage_per_amplicon                                1
_copy_mosdepth_summary_amplicons                                    1
_copy_mosdepth_timestep_coverage_per_amplicon_csv                   1
_copy_mosdepth_timestep_coverage_per_amplicon_png                   1
_copy_resources_yaml                                                1
_copy_samples_tsv                                                   1
_copy_sniffles2_var_calls                                           1
_copy_snv_indels_var_calls_ensembled_annotated_soft_filtered        1
_copy_soft_clipped_bam                                              1
_copy_soft_clipped_bam_bai                                          1
_copy_units_tsv                                                     1
_copy_yield_pool_1                                                  1
_copy_yield_pool_2                                                  1
_copy_yield_pool_3                                                  1
aligning_bam_softclip                                               1
aligning_bam_sort                                                   1
all                                                                 1
bam2fastq                                                           1
cnv_sv_bgzip                                                        1
cnv_sv_sniffles2_call                                               1
cnv_sv_tabix                                                        2
compress_fastq                                                      1
copy_annotation_vep                                                 1
dorado_align                                                        1
filtering_bcftools_include_region                                   1
filtering_filter_vcf                                                1
filtlong                                                            1
mosdepth_merge                                                      1
mosdepth_merge_timestep                                             1
mosdepth_overlap                                                   25
mosdepth_overlap_timestep                                           1
plot_yield_timestep                                                 1
qc_mosdepth_amplicons                                               1
qc_mosdepth_bed_per_exon                                            1
rename_vaf_to_af                                                    1
results_report_xlsx                                                 1
snv_indels_bgzip                                                    2
snv_indels_tabix                                                    3
yield_per_pool                                                      3
total                                                              73

Select jobs to execute...
```

**Note about the rule `dorado_align`**

The reference genome hg38 is too large to be included in this repository. 
Therefore, a dummy small reference genome is used, which consists of a reference sequence for the chromosome chr17 only 
(where TP53 is located).
Consequently, while the alignment is actually performed, 
the output aligned BAM file is not sensible as it contains unmapped reads that would normally map to other chromosomes.
For that reason, the small integration test overwrites the aligned BAM file with a more sensible file with exactly the same reads
and that was preliminarily aligned against the reference genome hg38:
`.tests/integration/test_data/preprocessed/D25-test007_T_reads.ont_adapt_trim.filtered.aligned.bam`

**Note about omitted rules**

 * the rules pertaining to basecalling with Dorado are not part of this integration test for the following reasons:
   * the POD5 files are too large to be included in the repository;
   * basecalling is time-consuming and requires GPU resources that are not always available on the CI runners or on local OS.
   Instead, a subsampled basecalled BAM file is provided as input to the test pipeline in 
   `.tests/integration/basecalling/dorado_duplex/D25-test007_T_reads.ont_adapt_trim.bam`;
 * the rule `pycoqc` is not included in this integration test as it requires a sequencing summary file that is not provided here.
   Indeed, the sequencing summary file that was actually generated for the test sample is too large to be included in the repository,
   and it cannot be downsampled easily;
 * the rule `sequali` is not included in this integration test as it requires the output files from `pycoqc`. 
   The corresponding lines in `.tests/integration/config/output_files.yaml` are therefore commented out;
 * the rule `multiqc` is not included in this integration test as it requires the output files from `pycoqc` and from `sequali`.
   The corresponding lines in `.tests/integration/config/output_files.yaml` are therefore commented out;;
 * the rule `annotation_vep` is skipped. 
   For one thing, VEP is best run with large cache files from Ensembl that are not proper to be stored in a GitHub repository.
   Additionally, testing the rule `annotation_vep` was [already done](https://github.com/hydra-genetics/annotation/actions/runs/17585961582/job/49953876925) in the hydra-genetics module `annotation`.



## :rocket: Usage

To execute this pipeline on your own data, you must first [create the sample files](https://hydra-genetics.readthedocs.io/en/latest/run_pipeline/create_sample_files/) `samples.tsv` and `units.tsv`.

To run the pipeline, you can then use a command similar to the one used for the small integration test above.

Those commands and other preliminary steps to the execution can be gathered in a bash script.
Examples of such script for a HPC cluster with Slurm installed on it are provided in
`workflow/scripts/start_marvin_multisample.sh` (sequencing run on a MinION flowcell and multiplexed barcoded samples) and 
`workflow/scripts/start_marvin_single.sh` (sequencing run on a Flongle flowcell and a single non barcoded sample).

**NB**: multiplexed barcoded samples require a [user-defined sample sheet](https://software-docs.nanoporetech.com/dorado/latest/barcoding/sample_sheet/) to be basecalled and demultiplexed.
**NBB**: Flongle flowcells are to be discontinued by ONT, therefore the script `start_marvin_single.sh` is provided for legacy reasons only.

```bash
# create input files for hydra-genetics
# 
```



If need to remove `__pycache__` for some user-written package: `rm -rf workflow/scripts/futils/__pycache__`.

### Output files

Without the `--notemp` option, temporary files are removed after execution of the pipeline.
Log files and benchmark files for each rule are left in the directories.
The files that must remain after execution are listed in `config/output_files.yaml`,
you may modify this file to add or remove some output files.
The main output files in `Results` are summarized in the table below:

| Folder or file    | Description                                                                                             |
|-------------------|---------------------------------------------------------------------------------------------------------|
| `Results`         | Folder with results of the analysis as well as the configurations of the pipeline and the sample files. |
| `Results/Bed`     | Design files used for coverage computation and variant calling.                                         |
| `Results/Data`    | Aligned and soft-clipped BAM files, VCF files with indels and SVs.                                      |
| `Results/Reports` | Excel report with quality metrics and prefiltered variant calls.                                        |

Note about mosdepth_bed: using "amplicons.bed" calculates cov per chrom --> values for TP53 are much lower since many amplicons there

### Program versions

Default container: `docker://hydragenetics/common:3.0.0`

| Program     | Version    | Container                                                   | Note                                                                                                                                                                                                                                             |
|-------------|------------|-------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| dorado      | 1.1.1      | `docker://ontresearch/dorado:latest`                        | Used in multiple steps: basecalling, demultiplexing, trimming, alignment. The basecalling model used in Dorado must be downloaded separately. Singularity> dorado --version [2025-11-28 09:25:46.359] [info] Running: "--version" 1.1.1+e72f1492 |
| samtools    | 1.21       | `docker://hydragenetics/samtools:1.21`                      | -                                                                                                                                                                                                                                                |
| Filtlong    | 0.2.1      | `docker://quay.io/biocontainers/filtlong:0.2.1--hdcf5f25_4` | The `--min_length` and `--max_length` command line options must be available.                                                                                                                                                                    |
| Mosdepth    | 0.3.6      | `docker://hydragenetics/mosdepth:0.3.6`                     | -                                                                                                                                                                                                                                                |
| ClairS-TO   | 0.3.1      | `docker://hkubal/clairs-to:latest`                          | No scientific publication related to ClairS-TO has been found yet.                                                                                                                                                                               |
| DeepSomatic | 1.8.0      | `docker://google/deepsomatic:1.8.0`                         | The snakemake command that executes variant calling in the suitable Singularity container must have the binding argument `/usr/lib/locale/:/usr/lib/locale/`.                                                                                    |
| VarDict     | 1.8.3      | `docker://hydragenetics/vardict:1.8.3`                      | Originally designed for short-read data but has an "amplicon" mode that seems compatible with ONT data, many false positive for short indels with low VAF though.                                                                                |
| Sniffles2   | 2.0.7      | `docker://hydragenetics/sniffles2:2.0.7`                    | -                                                                                                                                                                                                                                                |
| VT          | 2015.11.10 | `docker://hydragenetics/vt:2015.11.10`                      | -                                                                                                                                                                                                                                                |
| VEP         | 111.0      | `docker://hydragenetics/vep:111.0`                          | -                                                                                                                                                                                                                                                |
| Picard      | 2.25.4     | `docker://hydragenetics/picard:2.25.4`                     | -                                                                                                                                                                                                                                                |

## :judge: Rule Graph

![rule_graph_reference](images/neville_mx_amplicon.png)

## References

Karst, S.M., Ziels, R.M., Kirkegaard, R.H. et al. 
High-accuracy long-read amplicon sequences using unique molecular identifiers with Nanopore or PacBio sequencing. 
Nat Methods 18, 165–169 (2021). https://doi.org/10.1038/s41592-020-01041-y

Whitford, W., Hawkins, V., Moodley, K.S. et al. 
Proof of concept for multiplex amplicon sequencing for mutation identification using the MinION nanopore sequencer. 
Sci Rep 12, 8572 (2022). https://doi.org/10.1038/s41598-022-12613-7

J. Park et al.
DeepSomatic: Accurate somatic small variant discovery for multiple sequencing technologies.
bioRxiv 2024.08.16.608331; doi: https://doi.org/10.1101/2024.08.16.608331 

Lai Z, Markovets A, Ahdesmaki M, Chapman B, Hofmann O, McEwen R, Johnson J, Dougherty B, Barrett JC, and Dry JR. 
VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research. 
Nucleic Acids Res. Volume 44, Issue 11, 20 June 2016, Page e108.
doi: https://doi.org/10.1093/nar/gkw227

Oxford Nanopore PLC. 2023. Dorado. [GitHub repository](https://github.com/nanoporetech/dorado). 
Retrieved from https://github.com/nanoporetech/dorado