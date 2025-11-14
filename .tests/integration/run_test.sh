#!/bin/bash
# This script is used to run integration tests for the project.
# Usage: cd .tests/integration && bash ./run_test.sh

TEST_SAMPLE="D25-test007_T"

# Set up the folder structure as if basecalling would have been run
mkdir -p basecalling/dorado_duplex && cp test_data/bam_pass/reads.basecalled.bam basecalling/dorado_duplex/${TEST_SAMPLE}_reads.basecalled.bam

set -e
snakemake -s ../../workflow/Snakefile -n --configfile config/config.yaml --config runfolder=.tests/integration/test_data --use-singularity --singularity-args  "--cleanenv"

# Draw DAG and rulegraph --> workflow/scripts
# snakemake -s workflow/Snakefile --configfile config/config.yaml --config runfolder=.tests/integration/test_data --dag | dot -Tpdf > images/dag.pdf
# snakemake -s workflow/Snakefile --configfile config/config.yaml --config runfolder=.tests/integration/test_data --rulegraph | dot -Tpdf > images/rulegraph.pdf