runFolder=/projects/wp4/nobackup/ONT_dev_projects/CGU_2024_05_Amplicons_Hemato/CGU_2024_05_M24_2_251017 batchId=M24_2 flowcellId=FBD39643 runId=20251017_1440_MN48987_FBD39643_37b961fb sampleSheet=/projects/wp4/nobackup/ONT_dev_projects/CGU_2024_05_Amplicons_Hemato/CGU_2024_05_M24_2_251017/M24_2/20251017_1440_MN48987_FBD39643_37b961fb/Samplesheet_M24_2.csv

# Start pipeline
mkdir -p images
snakemake --profile profiles/slurm/ -s workflow/Snakefile --configfile config/config.yaml --config runfolder=${runFolder} batchid=${batchId} runid=${runId} multisample=True samplesheet="${sampleSheet}" --forceall --rulegraph > images/neville_mx_amplicon.dot && dot -Tpdf images/neville_mx_amplicon.dot > images/neville_mx_amplicon.pdf && dot -Tpng images/neville_mx_amplicon.dot > images/neville_mx_amplicon.png
