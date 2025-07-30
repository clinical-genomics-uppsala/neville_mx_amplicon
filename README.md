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

The module consists of ... build upon [hydra-genetics]().
Multiplexed amplicon sequencing.

Aim with the multiplexed setup: decrease the workload of the amplification step in the wet lab.
Caution: if two amplicons cover regions that overlap, those amplicons must be assigned to different pools.

The data were obtained from sequencing the prepared library on a MinION machine with the sequencing kit 14.
The duplex basecalling is performed with [Dorado](https://github.com/nanoporetech/dorado) 
([documentation](https://dorado-docs.readthedocs.io/en/latest/)).
Dorado is the currently recommended software for MinION output.

### Steps performed in the analysis

The main processing steps are:

- Basecalling with *Dorado* from POD5 files (raw output of ONT sequencing machines)
- Filtering with *Filtlong*: 
only reads with sufficient quality (based on the Phred quality scores) are kept. No reference is used. 
Reads that are too short or too long are filtered out.
Very long reads (over 4kb in our case) can be chimeras that result from concatenated amplicons.
The `--split` option could be used as alternative strategy, in such case a reference genome is required.
The thresholds for minimal and maximal length may be adjusted depending on the set of amplicons.
See the [documentation for Filtlong](https://github.com/rrwick/Filtlong/tree/main) for more options.
- Alignment with *Dorado* and soft-clipping with *samtools*,
- Variant calling with:
  - the software [ClairS-TO](https://github.com/HKU-BAL/ClairS-TO) via the Docker container provided by the development team,
  - the software [DeepSomatic](https://github.com/google/deepsomatic) via Docker container provided by the development team and the [documented examples](https://github.com/google/deepsomatic/blob/r1.8/docs/deepsomatic-case-study-ont-tumor-only.md),
  - the software VarDict [](),
  - the software Sniffles v2 [](),
- Annotation of the variants with VEP,
- Computation of some quality metrics:
  - Read quality with Sequali,
  - Estimated number of reads for each amplicon: approximated by the mean coverage computed with *mosdepth*,
  - Counts and proportion of reads for each amplicon relatively to the pool it belongs to,

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

### Set up the virtual environment

Create a virtual environment and install the required Python packages:
```bash
cd <root_of_the_repository>
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```


## :white_check_mark: Testing (TODO)

### Linting
Run the linter to check the code style and the correctness of the Snakefile:
```bash
cd <root_of_the_repository>
snakemake -s workflow/Snakefile -n --configfile .tests/integration/config/config.yaml --config runfolder=.tests/integration/test_data --debug-dag --use-singularity --singularity-args  " --cleanenv"
```

### Dry run
Configure the virtual environment `.venv` and test the pipeline with a dry run to see if the rules are correctly 
defined and the dependencies are satisfied:

```bash
cd <root_of_the_repository>
source .venv/bin/activate
snakemake -s workflow/Snakefile -n --configfile .tests/integration/config/config.yaml --config runfolder=.tests/integration/test_data --debug-dag --use-singularity --singularity-args  " --cleanenv"
```

NB: Add the option `--profile profiles/slurm/` to run the pipeline on a cluster with SLURM installed on it.

### Small integration test
Also serves as example of how to run the pipeline with a small test dataset.
The workflow repository contains a small test dataset `.tests/integration` which can be run like ...

```bash
...
```

Note on how the test data was generated:
```bash
$ cd .tests/integration/basecalling/dorado_duplex
$ samtools split D25-test007_T_reads.basecalled.bam
$ mv D25-test007_T_reads.basecalled_0.bam ../../test_data/bam_pass/ABC123_pass_456_789_0.bam
$ mv D25-test007_T_reads.basecalled_1.bam ../../test_data/bam_pass/ABC123_pass_456_789_1.bam
$ samtools index ../../test_data/bam_pass/ABC123_pass_456_789_0.bam
$ samtools index ../../test_data/bam_pass/ABC123_pass_456_789_1.bam
```
The file name "ABC123_pass_456_789" was chosen to match ONT's naming for the basecalled files.
BAI index files are usually part of the output data from the sequencer, but here we create them manually.


## :rocket: Usage

To use this module in your workflow, follow the description in the
[snakemake docs](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules).
Add the module to your `Snakefile` like so:

```bash
module prealignment:
    snakefile:
        github(
            "pipeline_pool_amplicon",
            path="workflow/Snakefile",
            tag="1.0.0",
        )
    config:
        config


use rule * from tmp.pipeline_pool_amplicon as pipeline_pool_amplicon_*
```

If need to remove `__pycache__` for some user-written package: `rm -rf workflow/scripts/futils/__pycache__`.

### Output files

The following output files should be targeted via another rule:

| File | Description |
|---|---|
| `pipeline_pool_amplicon/PATH/FILE` | DESCRIPTION |

### Program versions

Default container: ...

| Program     | Version | Container                                                 | Note                                                                                                                                                         |
|-------------|---------|-----------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------------|
| dorado      | 0.7.2   | `docker://ontresearch/dorado:latest`                      | Server version?                                                                                                                                              |
| samtools    | 1.21    | `docker://hydragenetics/samtools:1.21`                    | -                                                                                                                                                            |
| Filtlong    | 0.2.1   | `docker://quay.io/biocontainers/filtlong:0.2.1--hdcf5f25_4` | The `--min_length` and `--max_length` command line options must be available.                                                                                |
| Mosdepth    | 0.3.6   | `docker://hydragenetics/mosdepth:0.3.6`                   | -                                                                                                                                                            |
| ClairS-TO   | 0.3.1   | `docker://hkubal/clairs-to:latest`                     | No scientific publication related to ClairS-TO has been found yet.                                                                                           |
| DeepSomatic | 1.8.0   | `docker://google/deepsomatic:1.8.0`                   | The snakemake command that executes variant calling in the suitable Singularity container must have the binding argument `/usr/lib/locale/:/usr/lib/locale/` |


## :judge: Rule Graph

![rule_graph_reference](images/rulegraph.pdf)

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

Lai Z, Markovets A, Ahdesmaki M, Chapman B, Hofmann O, McEwen R, Johnson J, Dougherty B, Barrett JC, and Dry JR. VarDict: a novel and versatile variant caller for next-generation sequencing in cancer research. Nucleic Acids Res. 2016, pii: gkw227.

VarDict: https://academic.oup.com/nar/article/44/11/e108/2468301