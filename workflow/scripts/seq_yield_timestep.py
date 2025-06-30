import os
import sys
import subprocess
import pysam
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

from futils import *
import re

"""
Estimate the number of reads for each amplicon at increasing sequencing elapsed times.
Use mean coverage metrics from Mosdepth in targeted region to estimate the number of reads.
"""

indir = snakemake.input.indir
# indir = "results/mosdepth/timestep_coverage"
# indir = os.path.join(snakemake.config["proj_data"],
#                      snakemake.config["run_date"] + "_" + snakemake.config["run_id"],
#                      "bam_pass")

target_set = ["FLT3_ITD_3kb_1", "FLT3_TKD_3kb_1",
              "IDH1_1", "IDH2_1",
              "NPM1_5",
              "TP53_3kb_A1_only", "TP53_3kb_B1_only", "TP53_C10_only", "TP53_D2+J3", "TP53_D2_only",
              "TP53_E1_only", "TP53_F2_only", "TP53_G4_only", "TP53_H4_only", "TP53_I4_only"
              ]
# pool P2 only
p2 = ["IDH1_1", "TP53_3kb_B1_only", "TP53_G4_only", "TP53_D2+J3", "NPM1_5"]

batches = range(145)

readcounts = []
prev_counts = pd.Series(0.0, index=target_set, name="prev_counts")
prev_counts.index.name = "target"
for i, batch in enumerate(batches):
    dfcov = pd.read_csv(os.path.join(indir, f"timestep{batch}_coverage_per_amplicon.csv")).set_index("target")
    dfcov = dfcov.join(prev_counts, how="outer").drop("prev_counts", axis=1)  # expand axis in case an amplicon has not produced any read in the time interval
    dfcov = dfcov.loc[target_set].fillna(0.0)
    dfcov["timestep"] = batch * 10
    dfstep = dfcov[["timestep", "mean"]]
    readcounts.append(dfcov)

df = pd.concat(readcounts).reset_index(drop=False).set_index(["target", "timestep"])# [target_set]
dfcum = df.groupby(level=0).cumsum().reset_index(drop=False)  # cumulative sum over time per amplicon
dfcum.to_csv(os.path.join(snakemake.output.csv), index=False)
# dfcum.to_csv(os.path.join(indir, "cumsum_coverage_per_amplicon.csv"), index=False)

sns.set_style("whitegrid")
fig, ax = plt.subplots(1, 1, figsize=(16, 8))
sns.lineplot(data=dfcum, x="timestep", y="mean",
             hue="target", hue_order=target_set, palette="tab10", style="target",  # hue_order=target_set
             ax=ax)
ax.set_title("Estimated reads counts per amplicon")
ax.set_ylabel("Estimated reads count from mean coverage")
ax.set_xlabel("Sequencing time elapsed (in minutes)")
plt.savefig(os.path.join(snakemake.output.png))
# plt.savefig(os.path.join(indir, "cumsum_coverage_per_amplicon.png"))
plt.show()
