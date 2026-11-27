# Small integration test

The small integration test shows how to run the pipeline and how the working directory should look.
It runs on a laptop in about 2 minutes (tested on 22 cores, 32 GiB RAM, Ubuntu 22.04).

Some container images pulled during the test are larger than 8 GB, which exceeds the disk space of GitHub-hosted runners,
so the test is not run automatically in CI. Run it locally as described below.

## Test data

Sample `D25-test007` comes from a commercial [myeloid gDNA reference standard](https://horizondiscovery.com/en/reference-standards/products/myeloid-gdna-reference-standard)
(Horizon Discovery) with known variants in *IDH1*, *IDH2*, *NPM1*, *FLT3* and *TP53*.

!!! warning
    This reference standard is not optimised for ONT sequencing, and the results are not clinically meaningful.
    The test only checks that the pipeline runs end to end.

`.tests/integration/` contains:

* `test_data/bedfiles/`: the BED files of the design;
* `test_data/Mtest/ABC123/`: a run folder in MinKNOW layout, with BAM files in `bam_pass/` and an empty `pod5/`;
* `test_data/preprocessed/`: pre-computed files that replace steps that cannot run in the test (see below);
* `reference/`: a small reference made of chromosome 17 only (`TP53_chr17_GRCh38.fasta.gz`);
* `config/`: test configuration, which overrides `config/config.yaml`.

## Steps replaced or skipped

| Step | What happens in the test |
|------|--------------------------|
| Basecalling | Skipped: POD5 files are too large and need a GPU. A basecalled BAM file with about 15,000 reads is used as input. |
| Alignment | Runs against the chr17-only reference, then the output is replaced with a BAM file containing the same reads aligned to the whole of hg38 (`runid: ABC123` triggers this). |
| VEP annotation | Replaced by `copy_annotation_vep`, which copies pre-annotated VCF files from `test_data/preprocessed/`. |
| pycoQC, Sequali, MultiQC | Not run: the sequencing summary is too large to include. The matching lines in `.tests/integration/config/output_files.yaml` are commented out. |

## Preparation

```bash
# from the repository root, with the virtual environment set up
cp .tests/integration/config/*.tsv .tests/integration/     # if samples.tsv/units.tsv are missing
mkdir -p .tests/integration/basecalling/dorado_duplex/
cp .tests/integration/test_data/preprocessed/Mtest_D25-test007_T_reads.ont_adapt_trim.bam \
   .tests/integration/basecalling/dorado_duplex/
mkdir -p .tests/integration/tmp/
bgzip -d .tests/integration/reference/TP53_chr17_GRCh38.fasta.gz
```

The machine needs internet access to pull the container images.

## Execution

```bash
cd .tests/integration/
source ../../.venv/bin/activate
snakemake -s ../../workflow/Snakefile -j 1 --show-failed-logs \
    --configfiles ../../config/config.yaml config/config.yaml \
    --use-singularity --singularity-prefix .snakemake/singularity/ \
    --singularity-args " --cleanenv --containall --bind $PWD/tmp:/tmp -B $PWD:$PWD -B $HOME -B /usr/lib/locale/:/usr/lib/locale/ --disable-cache "
```

Adapt `--singularity-args` to your system if needed. Add `--notemp` to keep the intermediate files.

## Results

The main output is written to `.tests/integration/Results/` (see [Result files](result_files.md)).

## Other checks

Lint the workflow:

```bash
cd .tests/integration
snakemake --lint -n -s ../../workflow/Snakefile --configfiles config/config.yaml
```

Dry run:

```bash
cd .tests/integration
snakemake -n -s ../../workflow/Snakefile --configfile config/config.yaml --config runfolder=../../.tests/integration/test_data
```

The GitHub Actions workflows in `.github/workflows/` run linting, `snakefmt`, `pycodestyle`, `pytest` and a dry run on each pull request.
