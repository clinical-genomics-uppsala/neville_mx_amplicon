#!/usr/bin/env bash
# To run script: ...

# set -euox pipefail

module load slurm-drmaa/1.1.4
module load singularity/3.7.1
module load snakemake/7.22.0
module load samtools/1.17

# Parse command line arguments
projFolder=/beegfs-storage/projects/wp4/nobackup/workspace/camille_test/ampliconthemato/neville_mx_amplicon
runFolder=/projects/wp4/nobackup/ONT_dev_projects/CGU_2024_05_Amplicons_Hemato/CGU_2024_05_MWash1_250804
batchId=MWash1
# runFolder=/projects/wp4/nobackup/ONT_dev_projects/CGU_2024_05_Amplicons_Hemato/CGU_2024_05_M21_wash2_250806
# batchId=M21
# runFolder=/projects/wp4/nobackup/ONT_dev_projects/CGU_2024_05_Amplicons_Hemato/CGU_2024_05_Mwash3_250811
# batchId=Mwash3
# runFolder=/projects/wp4/nobackup/ONT_dev_projects/CGU_2024_05_Amplicons_Hemato/CGU_2024_05_Wash4_250813
# batchId=Wash4
flowcellId=FBB06783

runId=$( ls -1 "$runFolder/${batchId}" | grep ${flowcellId} )
echo "runfolder: '${runFolder}/${batchId}/${runId}'" > runfolder.txt
sampleSheet=${runFolder}/${batchId}/${runId}/SAMPLESHEET_ONT_MWASH.csv
csvDelim=','


# Merge BAM files and p per sample and create input files for the pipeline
while IFS=$csvDelim read -r position_id flow_cell_id kit experiment_id sample_id alias barcode; do
  echo "$sample_id has barcode $barcode."
  mkdir -p ${runFolder}/${sample_id}/${runId}/bam_pass_merged
  echo "Merging BAM files found for $sample_id into ${runFolder}/${sample_id}/${runId}/bam_pass_merged"
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
  cd ${projFolder}
  # Set merged BAM file as input for the pipeline to start
  mkdir -p basecalling/dorado_duplex
  cp ${runFolder}/${sample_id}/${runId}/bam_pass_merged/reads.basecalled.bam ./basecalling/dorado_duplex/${sample_id}_T_reads.basecalled.bam
  # Create input files to the pipeline per sample
  source .venv/bin/activate
  hydra-genetics create-input-files -d ${runFolder}/${sample_id}/${runId}/bam_pass_merged/ -t T -p ONT -f
  # --default-barcode $barcode
  cp units.tsv units_$alias.tsv
  cat samples.tsv | cut -d$'\t' -f1 > samples_$alias.tsv
done < <(tail -n +2 ${sampleSheet}) # skip header line while reading csv

cat units.tsv | head -1 > header_units.tsv
cat samples.tsv | cut -d$'\t' -f1 | head -1 > header_samples.tsv

cd ${projFolder}
mkdir -p tmp
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

# Start pipeline
snakemake --profile profiles/slurm/ -s workflow/Snakefile \
--configfile config/config.yaml \
--config runfolder=${runFolder} batchid=${batchId} runid=${runId} multisample=True samplesheet=${runFolder}/${batchId}/${runId}/SAMPLESHEET_ONT_MWASH.csv \
--notemp
