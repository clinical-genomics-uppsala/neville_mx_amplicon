#!/usr/bin/env bash

# Example of usage: $ bash workflow/scripts/fastq_to_rehead_bam.sh /projects/wp4/nobackup/ONT_dev_projects/CGU_2024_05_Amplicons_Hemato/CGU_2025_22_exp2_251212_RO CGU_2025_22_exp2 FBE72298 pass

module load samtools/1.17

# Parse command line arguments
projFolder=/beegfs-storage/projects/wp4/nobackup/workspace/camille_test/ampliconthemato/neville_mx_amplicon
runFolder=$1
batchId=$2
flowcellId=$3
filePassFail=$4

# runFolder=/projects/wp4/nobackup/ONT_dev_projects/CGU_2024_05_Amplicons_Hemato/CGU_2025_22_exp2_251212_RO
# batchId=CGU_2025_22_exp2
# flowcellId=FBE72298

runId=$( ls -1 "$runFolder/${batchId}" | grep ${flowcellId} )
sampleSheet=${runFolder}/${batchId}/${runId}/Samplesheet_${batchId}.csv
csvDelim=','

# Convert fastq pass/fail files to uBAM, reheader the uBAM, merge time-stepped BAM files per sample, and create input files for the pipeline
while IFS=$csvDelim read -r position_id flow_cell_id kit experiment_id sample_id alias barcode; do
# strip whitespace characters at the end of the CSV fields
  barcode="${barcode%[[:space:]]}"
  # save sample_id to avoid overwriting it in the for-loop downstreams
  sampleId=${sample_id}
  echo "${sampleId} has barcode $barcode."
  mkdir -p ${runFolder}/${sampleId}/${runId}/bam_${filePassFail}
  fastq_files=$( ls ${runFolder}/${batchId}/${runId}/fastq_${filePassFail}/${barcode}/*.fastq.gz )
  echo "Converting FASTQ files to uBAM files for ${sampleId} into ${runFolder}/${sampleId}/${runId}/bam_${filePassFail}"
  
  for fastq in $fastq_files; do
    echo "Processing ${fastq}"
    basename=$( basename ${fastq} | cut -d'.' -f1 )
    METAINFO=$( zcat ${fastq} | head -1 | cut -d ' ' -f2-12 )
    eval $METAINFO
    # /!\ Now ${sample_id} = ${batchId}
    BAM_HEADER='@RG\tID:'${runid}'_'${basecall_model_version_id}'_'${barcode}'\tDT:'${start_time}'\tDS:runid='${runid}' basecall_model='${basecall_model_version_id}'\tLB:'${sample_id}'\tPL:ONT\tPM:'${position_id}'\tPU:'${flow_cell_id}'\tSM:'${barcode}'\tal:'${barcode_alias}
    samtools import --rg-line "$(echo -e ${BAM_HEADER})" $fastq -o ${runFolder}/${sampleId}/${runId}/bam_${filePassFail}/${basename}.bam
    samtools index ${runFolder}/${sampleId}/${runId}/bam_${filePassFail}/${basename}.bam
  done
done < <(tail -n +2 ${sampleSheet}) # skip header line while reading csv
