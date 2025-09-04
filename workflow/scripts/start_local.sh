#!/usr/bin/env bash
# To run script:
# bash start_marvin_single.sh... <path_to_the_runFolder>

set -euo pipefail

projFolder=/home/camille/ampliconthemato/pipeline_pool_amplicon
# runFolder=/projects/wp4/nobackup/ONT_dev_projects/CGU_2024_05_Amplicons_Hemato/CGU_2024_05_PoolM12_M13_RO/D21-03300/20250117_1351_MN48987_AWY831_f6f7634a
# runFolder=/projects/wp4/nobackup/ONT_dev_projects/CGU_2024_05_Amplicons_Hemato/CGU_2024_05_PoolM12_M13_RO
# runid=$( ls $runFolder/D21-03300 | grep AWY831 )
# echo "runfolder: '$runFolder/D21-03300/$runid'" >> config/config.yaml
# runFolder=$1

## Prepare the uBAM data to be able to create input files to the pipeline
## Merged uBAM must be in a terminal subfolder (no child directories) because create-input-files recursively searches for BAM files
#do the following if no bam pass merged yet and encapsulate in sbatch
# mkdir ${runFolder}/bam_pass_merged
# samtools merge ${runFolder}/bam_pass_merged/reads.basecalled.bam ${runFolder}/bam_pass/*.bam

# Start pipeline
cd ${projFolder}
source .venv/bin/activate
# hydra-genetics create-input-files -d ${runFolder}/bam_pass_merged/ -t T -p ONT -f
snakemake --cores 16 --profile profiles/local/ -s workflow/rules/dummy.smk \
--configfile config/config.yaml --notemp
# --dry-run
