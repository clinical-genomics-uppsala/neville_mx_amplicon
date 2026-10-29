#!/bin/bash
# Align the unfiltered reads of one sample and soft-clip the primers, mirroring the pipeline's
# alignment_ont_dorado_align -> alignment_ont_bam_sort/index -> alignment_ont_bam_softclip for the filtered reads.
# The output BAM is the --unfiltered-bam input of tools/noise_viz.
#
# Usage: sbatch align_unfiltered.sbatch /path/to/<sample>.ont_reads.unfiltered.fastq.gz [output prefix]
#
#SBATCH --job-name=align_unfiltered
#SBATCH --nodelist=compute10
#SBATCH --partition=core
#SBATCH --cpus-per-task=20
#SBATCH --mem=120G
#SBATCH --time=02:00:00
#SBATCH --output=align_unfiltered_%j.log

set -euo pipefail

FQ=${1:?usage: sbatch align_unfiltered.sbatch <unfiltered.fastq.gz> [output prefix]}
OUT=${2:-$(basename "$FQ" .fastq.gz).aligned}

# Paths from the pipeline config used for M20_4 (ref_data, bed_files): check the primer BED matches the pool
REF=/data/ref_genomes/GRCh38/GCA_000001405.15/GCA_000001405.15_GRCh38_no_alt_analysis_set_masked_chr.fasta
PRIMERS=/beegfs-storage/projects/wp4/nobackup/workspace/camille_test/ampliconthemato/bedfiles/M5-6-7/amplicons-primers.bed

DORADO=docker://hydragenetics/dorado:1.1.1
SAMTOOLS=docker://hydragenetics/samtools:1.21
BIND="-B /data,/beegfs-storage,/projects,$PWD,$(dirname "$(realpath "$FQ")")"
THREADS=${SLURM_CPUS_PER_TASK:-20}

echo "$(date) align $FQ -> $OUT.bam"
apptainer exec $BIND $DORADO dorado aligner -t "$THREADS" "$REF" "$FQ" > "$OUT.bam"

echo "$(date) sort and index"
apptainer exec $BIND $SAMTOOLS samtools sort -@ "$THREADS" -o "$OUT.sorted.bam" "$OUT.bam"
apptainer exec $BIND $SAMTOOLS samtools index "$OUT.sorted.bam"

echo "$(date) soft-clip primers"
apptainer exec $BIND $SAMTOOLS samtools ampliconclip \
    -o "$OUT.sorted.soft-clipped.bam" -b "$PRIMERS" -f "$OUT.ampliconclip.txt" \
    --soft-clip --both-ends --clipped "$OUT.sorted.bam"
apptainer exec $BIND $SAMTOOLS samtools index "$OUT.sorted.soft-clipped.bam"

rm "$OUT.bam"
echo "$(date) done: $OUT.sorted.soft-clipped.bam"
