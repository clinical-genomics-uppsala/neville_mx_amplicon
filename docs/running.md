# Running the Pipeline

This section describes how to validate the pipeline rules and run the small integration test locally.

---

## 1. Code Linting

To check the pipeline code style and structure, activate the virtual environment and run the Snakemake linter from the integration test folder:
```bash
source .venv/bin/activate
cd .tests/integration
snakemake --lint -n -s ../../workflow/Snakefile --configfiles config/config.yaml
```

---

## 2. Dry Run

A dry run validates that the rules, input files, configurations, and outputs are correctly defined without launching actual computations:
```bash
source .venv/bin/activate
cd .tests/integration
snakemake -n -s ../../workflow/Snakefile --configfile config/config.yaml --config runfolder=../../.tests/integration/test_data
```

---

## 3. Running the Integration Test

The integration test operates on a downsampled test sample (`D25-test007`) containing known variants in `IDH1`, `IDH2`, `NPM1`, `FLT3`, and `TP53`.

### Preparations
1. Create the `tmp/` folder within the integration test directory:
   ```bash
   mkdir -p .tests/integration/tmp/
   ```
2. Decompress the reference genome if necessary:
   ```bash
   bgzip -d .tests/integration/reference/TP53_chr17_GRCh38.fasta.gz
   ```
3. Set up the input test BAM file:
   ```bash
   mkdir -p .tests/integration/basecalling/dorado_duplex/
   cp .tests/integration/test_data/preprocessed/Mtest_D25-test007_T_reads.ont_adapt_trim.bam .tests/integration/basecalling/dorado_duplex/
   ```

### Execution
Execute the pipeline within the container environment using Singularity:
```bash
cd .tests/integration/
source ../../.venv/bin/activate
snakemake -s ../../workflow/Snakefile -j 1 --show-failed-logs \
  --configfiles ../../config/config.yaml config/config.yaml \
  --use-singularity --singularity-prefix .snakemake/singularity/ \
  --singularity-args "--cleanenv --containall --bind $PWD/tmp:/tmp -B $PWD:$PWD -B $HOME -B /usr/lib/locale/:/usr/lib/locale/ --disable-cache"
```
*Note: The test takes approximately 2 minutes to complete on a modern CPU.*
