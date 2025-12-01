#!/usr/bin/env bash
# To run script:
# bash start_marvin_single.sh... <path_to_the_runFolder> <sample_ID> <flowcell_ID>
# Example: bash workflow/scripts/start_marvin_single.sh
# /projects/wp4/nobackup/ONT_dev_projects/CGU_2024_05_Amplicons_Hemato/CGU_2024_05_PoolM10 D21-06478 AWJ762

set -euox pipefail

module load slurm-drmaa/1.1.4
module load singularity/3.7.1
module load snakemake/7.22.0
module load samtools/1.17

projFolder=/beegfs-storage/projects/wp4/nobackup/workspace/camille_test/ampliconthemato/neville_mx_amplicon
runFolder=$1
sampleId=$2
flowcellId=$3
runId=$( ls -1 "$runFolder/${sampleId}" | grep ${flowcellId} )

## Prepare the uBAM data to be able to create input files to the pipeline
## Merged uBAM must be in a terminal subfolder (no child directories) because create-input-files recursively searches for BAM files
# do the following if no bam pass merged yet
mkdir -p "${runFolder}/${sampleId}/${runId}/bam_pass_merged"
if [ ! -f "${runFolder}/${sampleId}/${runId}/bam_pass_merged/reads.basecalled.bam" ]
then
  cd ${runFolder}/${sampleId}/${runId}/bam_pass/
  ls -1 . | grep -iE '.+bam$' > "bam_list.txt"
  samtools merge -o ../bam_pass_merged/reads.basecalled.bam -b bam_list.txt
fi
cd ${projFolder}

# Start pipeline
cd ${projFolder}
mkdir -p tmp
source .venv/bin/activate
hydra-genetics create-input-files -d ${runFolder}/${sampleId}/${runId}/bam_pass_merged/ -t T -p ONT -f

# Measure running time
STARTTIME=$(date +%s)

snakemake --profile profiles/slurm/ -s workflow/Snakefile \
--configfile config/config.yaml \
--config runfolder=${runFolder} batchid=${sampleId} runid=${runId} multisample=False \
--notemp

ENDTIME=$(date +%s)
RUNTIME=$((ENDTIME - STARTTIME))
echo "Pipeline finished in $RUNTIME seconds."