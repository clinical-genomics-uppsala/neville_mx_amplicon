import os
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

"""
Estimate the number of reads for each amplicon at increasing sequencing elapsed times.
Use mean coverage metrics from Mosdepth in targeted region to estimate the number of reads.
"""

indir = snakemake.input.indir

target_set = ["FLT3_ITD_3kb_1", "FLT3_TKD_3kb_1",
              "IDH1_1", "IDH2_1",
              "NPM1_5",
              "TP53_3kb_A1_only", "TP53_3kb_B1_only", "TP53_C10_only", "TP53_D2+J3", "TP53_D2_only",
              "TP53_E1_only", "TP53_F2_only", "TP53_G4_only", "TP53_H4_only", "TP53_I4_only"
              ]
target_set_j3 = target_set + ["TP53_J3_only"]  # add J3 amplicon to the target set
target_set_j3.remove("TP53_D2+J3")

# pool P2 only
p2 = ["IDH1_1", "TP53_3kb_B1_only", "TP53_G4_only", "TP53_D2+J3", "NPM1_5"]

batches = range(145)

readcounts = []
prev_counts = pd.Series(0.0, index=target_set, name="prev_counts")
prev_counts.index.name = "target"
for i, batch in enumerate(batches):
    try:
        dfcov = pd.read_csv(os.path.join(indir, f"timestep{batch}_coverage_per_amplicon.csv")).set_index("target")
        # expand axis in case an amplicon has not produced any read in the time interval
        dfcov = dfcov.join(prev_counts, how="outer").drop("prev_counts", axis=1)
        dfcov = dfcov.loc[target_set].fillna(0.0)
        dfcov["timestep"] = (batch + 1) * 60  # * 10 if Flongle and * 60 if MinION
        dfstep = dfcov[["timestep", "mean"]]
        readcounts.append(dfcov)
    except FileNotFoundError:  # if seqrun failed or is stopped before elapsed time of 24h, then there < 144 bam files
        continue

df = pd.concat(readcounts).reset_index(drop=False).set_index(["target", "timestep"])
dfcum = df.groupby(level=0).cumsum().reset_index(drop=False)  # cumulative sum over time per amplicon
dfcum.to_csv(os.path.join(snakemake.output.csv), index=False)
# dfcum.to_csv(os.path.join(indir, "cumsum_coverage_per_amplicon.csv"), index=False)
print(dfcum[dfcum["target"].isin(["TP53_D2+J3", "TP53_D2_only"])])
for t in set(sorted(dfcum["timestep"])):
    row_j3 = dict(zip(dfcum.columns, [
        "TP53_J3_only",
        t,
        0,
        0,
        dfcum[(dfcum["target"] == "TP53_D2+J3") & (dfcum["timestep"] == t)]["mean"].values[0]
        - dfcum[(dfcum["target"] == "TP53_D2_only") & (dfcum["timestep"] == t)]["mean"].values[0]
    ]))
    dfcum = pd.concat([dfcum,
                       pd.DataFrame(data=row_j3, index=[0])],
                      axis=0, ignore_index=True)
dfcum = dfcum[dfcum["target"] != "TP53_D2+J3"]  # remove the D2+J3 amplicon
print(dfcum)

sns.set_style("whitegrid")
fig, ax = plt.subplots(1, 1, figsize=(16, 8))
sns.lineplot(data=dfcum, x="timestep", y="mean",
             hue="target", hue_order=target_set_j3, palette="tab10", style="target",  # hue_order=target_set
             ax=ax)
ax.set_title("Estimated reads counts per amplicon")
ax.set_ylabel("Estimated reads count from mean coverage")
ax.set_xlabel("Sequencing time elapsed (in minutes)")
plt.savefig(os.path.join(snakemake.output.png))
# plt.savefig(os.path.join(indir, "cumsum_coverage_per_amplicon.png"))
plt.show()
