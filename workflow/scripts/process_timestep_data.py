import os
import subprocess

"""
This script creates and processes mosdepth summary files as well as it extracts relevant data for each target.
"""

for bamfile in os.listdir(snakemake.input.bamdir):
    print(f"\r\nProcessing BAM file: {snakemake.input.bamdir}/{bamfile}")
    if bamfile.endswith(".bam"):
        nbam = int(bamfile.split('_')[-1].replace(".bam", ""))
        for bedfile in snakemake.input.amplibed:
            target = os.path.basename(bedfile.replace(".bed", ""))
            nbamdir = os.path.join(snakemake.output.outdir, f"{nbam}")
            os.makedirs(nbamdir, exist_ok=True)
            print(f"Processing target: {target} at timestep {nbam}")
            mosdepth_cmd = (
                f"chrom=$( cat {bedfile} | cut -d$'\t' -f1 )"
                f" && "
                f"mosdepth"
                f" -t {snakemake.resources.threads}"
                f" -c $chrom"
                f" -b {bedfile}"
                f" {nbamdir}/{target}"
                f" {snakemake.input.bamdir}/{bamfile}"
                f" 2>> {snakemake.log}"
            )
            subprocess.run(mosdepth_cmd, shell=True, check=True, executable="/bin/bash")
