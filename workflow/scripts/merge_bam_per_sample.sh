#!/usr/bin/env bash
# To run script: ...

# set -euox pipefail

module load samtools/1.17

projFolder=/beegfs-storage/projects/wp4/nobackup/workspace/camille_test/ampliconthemato/neville_mx_amplicon
runFolder=/projects/wp4/nobackup/ONT_dev_projects/CGU_2024_05_Amplicons_Hemato/CGU_2024_05_MWash1_250804
batchId=MWash1
runId=20250804_1309_MN48987_FBB06783_fd3e24b1
sampleSheet=${runFolder}/${batchId}/${runId}/SAMPLESHEET_ONT_MWASH1.csv
csvDelim=$'\t'

while IFS=$csvDelim read -r position_id flow_cell_id kit experiment_id sample_id alias barcode; do
    echo "$sample_id has barcode $alias."
    mkdir -p ${runFolder}/${sample_id}/${runId}/bam_pass_merged
    echo "Merging BAM files found for $sample_id into ${runFolder}/${sample_id}/${runId}/bam_pass_merged"
    if [ ! -f "${runFolder}/${sample_id}/${runId}/bam_pass_merged/reads.basecalled.bam" ]
	then
	  cd ${runFolder}/${batchId}/${runId}/bam_pass/${alias}
	  ls -1 . | grep -iE '.+bam$' > "bam_list.txt"
	  samtools merge -o ${runFolder}/${sample_id}/${runId}/bam_pass_merged/reads.basecalled.bam -b bam_list.txt
	fi
	cd ${projFolder}
	source .venv/bin/activate
	hydra-genetics create-input-files -d ${runFolder}/${sample_id}/${runId}/bam_pass_merged/ -t T -p ONT -f
	# --default-barcode $barcode
	cp units.tsv units_$alias.tsv
	cp samples.tsv samples_$alias.tsv
done < <(tail -n +2 ${sampleSheet}) # skip header line while reading csv

cat units.tsv | head -1 > header_units.tsv
cat samples.tsv | head -1 > header_samples.tsv

cd ${projFolder}
mkdir -p tmp
# Merge tsv files together, assumes only 1 sample per file
cp header_units.tsv units.tsv
cp header_samples.tsv samples.tsv
while IFS=$csvDelim read -r position_id flow_cell_id kit experiment_id sample_id alias barcode; do
  echo -e $sample_id$'\t'$( cat units_$alias.tsv | tail -1 | cut --complement --output-delimiter='\t' -d$'\t' -f1 ) >> units.tsv
	echo -e $sample_id$'\t'$( cat samples_$alias.tsv | tail -1 | cut --complement --output-delimiter='\t' -d$'\t' -f1 ) >> samples.tsv
	rm -f units_$alias.tsv
  rm -f samples_$alias.tsv
done < <(tail -n +2 ${sampleSheet})

# Start pipeline
snakemake --profile profiles/slurm/ -s workflow/Snakefile \
--configfile config/config.yaml --config runfolder=${runFolder}/${sampleId}/${runId} --notemp