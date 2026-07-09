#!/usr/bin/env bash

# Usage: create directory for analysis: $ mkdir MXA999 && cd MXA999
# Usage: run this script: $ bash /proj/ngi2024001/nobackup/bin/wp2_mxa/start_wp2_mx_amplicon_minion.sh <runFolder> <batchId> <flowcellId> [OTHER_PARAMS] [...]

# Directory where the MinKNOW data are saved, the samplesheet should be there as well
runFolder=$1
# Short name/ID of the sequencing experiment
batchId=$2
# ID of the flowcell that was used for the sequencing run
flowcellId=$3

set -eox pipefail

echo "RUNNING: wp2 neville_mx_amplicon (MXA)"

# git_repo_url_neville="https://github.com/clinical-genomics-uppsala/neville_mx_amplicon.git"

# Initialize variables
# bin_path="/projects/bin/wp2_mxa"
bin_path="/proj/ngi2024001/nobackup/bin/wp2_mxa"
inbox_path=""
analysis_path=$(pwd)
sequenceid=""

projFolder=${analysis_path}
runId=$( ls -1 "$runFolder/${batchId}" | grep ${flowcellId} )
sampleSheet=${runFolder}/${batchId}/${runId}/Samplesheet_${batchId}.csv
csvDelim=','

# Only change repo version here
# neville_version="develop"
neville_version="offline-pack"


# Process options and arguments
# TODO

# Add blank line at the end of the samplesheet if need be
sed -i '$a\' ${sampleSheet}

# Check if required options are provided
#if [ -z "$inbox_path" ]; then
#    echo "Error: --inbox-path is required."
#    exit 1
#fi


# If no sequence-id is defined, set it to current folder name
if [ -z "$sequenceid" ]; then
    sequenceid=${PWD##*/}
fi

neville_path=${bin_path}/${neville_version}
source ${neville_path}/venv/bin/activate
echo "Source code for the pipeline: ${neville_path}"

# Load required modules (Miarka)
module load slurm-drmaa &&
	module load bioinfo-tools samtools/1.17

## Prepare the uBAM data to be able to create input files to the pipeline
## Merged uBAM must be in a terminal subfolder (no child directories) because create-input-files recursively searches for BAM files
# Merge BAM files and p per sample and create input files for the pipeline
while IFS=$csvDelim read -r position_id flow_cell_id kit experiment_id sample_id alias barcode; do
  # strip whitespace characters at the end of the CSV fields
  barcode="${barcode%[[:space:]]}"
  echo "${sample_id} has barcode $barcode."
  mkdir -p ${runFolder}/${sample_id}/${runId}/bam_pass_merged
  echo "Merging BAM files found for ${sample_id} into ${runFolder}/${sample_id}/${runId}/bam_pass_merged"
  if [ ! -f "${runFolder}/${sample_id}/${runId}/bam_pass_merged/reads.basecalled.bam" ]
  then
    cd ${runFolder}/${batchId}/${runId}/bam_pass/${barcode}
    ls -1 . | grep -iE '.+bam$' > "bam_list.txt"
    samtools merge -o ${runFolder}/${sample_id}/${runId}/bam_pass_merged/reads.basecalled.bam -b bam_list.txt
  fi
  # Restructure time-stepped BAM files to be per sample
  cd ${runFolder}/${sample_id}/${runId}
  mkdir -p bam_pass
  rsync -ruv ${runFolder}/${batchId}/${runId}/bam_pass/${barcode}/* ./bam_pass/
  # Prep
  cd ${analysis_path}
  # Set merged BAM file as input for the pipeline to start
  mkdir -p basecalling/dorado_duplex
  cp ${runFolder}/${sample_id}/${runId}/bam_pass_merged/reads.basecalled.bam ./basecalling/dorado_duplex/${sample_id}_T_reads.basecalled.bam
  # Create input files to the pipeline per sample
  source ${neville_path}/venv/bin/activate
  hydra-genetics create-input-files -d ${runFolder}/${sample_id}/${runId}/bam_pass_merged/ -t T -p ONT -f
  cp units.tsv units_$alias.tsv
  cat samples.tsv | cut -d$'\t' -f1 > samples_$alias.tsv
done < <(tail -n +2 ${sampleSheet}) # skip header line while reading csv

# Edit samples.tsv and units.tsv to handle all samples
cat units.tsv | head -1 > header_units.tsv
cat samples.tsv | cut -d$'\t' -f1 | head -1 > header_samples.tsv

cd ${analysis_path}
mkdir -p tmp
# Copy configs for common.smk to work
cp -r ${neville_path}/neville_mx_amplicon/config ${analysis_path}/

# Merge tsv files together, assumes only 1 sample per file
cp header_units.tsv units.tsv
cp header_samples.tsv samples.tsv
while IFS=$csvDelim read -r position_id flow_cell_id kit experiment_id sample_id alias barcode; do
  echo -e $sample_id$'\t'$( cat units_$alias.tsv | tail -1 | cut --complement --output-delimiter='\t' -d$'\t' -f1 ) >> units.tsv
  # echo -e $sample_id$'\t'$( cat samples_$alias.tsv | tail -1 | cut --complement --output-delimiter='\t' -d$'\t' -f1 ) >> samples.tsv
  echo -e $sample_id >> samples.tsv
  rm -f units_$alias.tsv
  rm -f samples_$alias.tsv
  rm -f header_units.tsv
  rm -f header_samples.tsv
done < <(tail -n +2 ${sampleSheet})

# Measure running time
STARTTIME=$(date +%s)

# Start pipeline
snakemake --profile ${neville_path}/neville_mx_amplicon/profiles/miarka_slurm -s ${neville_path}/neville_mx_amplicon/workflow/Snakefile \
	--configfile config/config_miarka.yaml \
	--config runfolder=${runFolder} batchid=${batchId} runid=${runId} multisample=True samplesheet="${sampleSheet}" TAG_OR_BRANCH=${neville_version} \
	--notemp
ENDTIME=$(date +%s)
RUNTIME=$((ENDTIME - STARTTIME))
echo "Pipeline finished in $RUNTIME seconds."
