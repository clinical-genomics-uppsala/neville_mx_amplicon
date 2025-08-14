import os
import subprocess

# This script creates and processes mosdepth summary files as well as it extracts relevant data for each target.

# 1. Compute coverage in non overlapping regions for each target at each timestep
# input:
#     bamdir = os.path.join(config["runfolder"], "{sample}", config["runid"], "bam_pass"),
#     amplibed = [
#         f"{config.get('bed_files')}/{target}.bed"
#         for target in config.get("amplicons") + config.get("extra_regions")
#     ],
# output:
#     outdir = temp(directory("results/mosdepth/timestep/{sample}")),

for bamfile in os.listdir(snakemake.input.bamdir):
    print("\r\nProcessing BAM file:", bamfile)
    if bamfile.endswith(".bam"):
        nbam = int(bamfile.split('_')[-1].replace(".bam", ""))
        for bedfile in os.listdir(snakemake.input.amplibed):
            target = bedfile.replace(".bed", "")
            # Create output directory for the current timestep
            nbamdir = os.path.join(snakemake.output.outdir, f"{nbam}")
            os.makedirs(nbamdir, exist_ok=True)
            print(f"Processing target: {target} at timestep {nbam}")

            # # Define output files
            # out_bed = os.path.join(outdir, f"{target}.regions.bed.gz")
            # out_csi = out_bed + ".csi"
            # out_global = os.path.join(outdir, f"{target}.mosdepth.global.dist.txt")
            # out_region = os.path.join(outdir, f"{target}.mosdepth.region.dist.txt")
            # out_summary = os.path.join(outdir, f"{target}.mosdepth.summary.txt")
            #
            # # Log file
            # log = os.path.join(outdir, f"{target}.log")

                # Run mosdepth command
            mosdepth_cmd = (
                f"chrom=$( cat {bedfile} | cut -d$'\t' -f1 )"
                f" && "
                f"mosdepth"
                f" -t {snakemake.resources.threads}"
                f" -c $chrom"
                f" -b {bedfile}"
                f" {nbamdir}/{target}"
                f" {snakemake.input.bam}"
                f" 2>> {snakemake.log}"
            )
            subprocess.run(mosdepth_cmd, shell=True, check=True, executable="/bin/bash")

