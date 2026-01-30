import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import markers
from matplotlib import colorbar
from matplotlib.lines import Line2D
from matplotlib.font_manager import findfont, FontProperties
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.metrics import r2_score
import os

fig_format = "png"
vaf_csv = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/VAF_TM_ONT_260126.csv"
exp_samples_csv = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/list_exp_id_sample_id.csv"
cumcounts = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/cum_read_counts"

relabels = {
    'TP53_3kb_A1_only': 'TP53_A',
    'TP53_H4_only': 'TP53_H',
    'NPM1_5': 'NPM1',
    'TP53_3kb_B1_only': 'TP53_B',
    'TP53_G4_only': 'TP53_G',
    'TP53_F2_only': 'TP53_F',
    'TP53_I4_only': 'TP53_I',
    'TP53_C10_only': 'TP53_C',
    'TP53_E1_only': 'TP53_E',
    'FLT3_ITD_3kb_1': 'FLT3_ITD',
    'IDH2_1': 'IDH2',
    'TP53_J3': 'TP53_J',
    'IDH1_1': 'IDH1',
    'FLT3_TKD_3kb_1': 'FLT3_TKD',
    'TP53_D2_only': 'TP53_D',
}

df_samples = pd.read_csv(vaf_csv, sep=",")
df_exp_samples = pd.read_csv(exp_samples_csv, sep=",")
df_samples = (df_samples
              .merge(df_exp_samples, left_on="sample", right_on="sample_id", how="left")
              .drop(columns=["sample_id", "data_dir"]))

# Read cumulative counts after 1 hour

cumcounts_per_sample = []
for csv in os.listdir(cumcounts):
    if csv[-4:] == ".csv" and csv.split("_")[1] == "MinION":
        dfrun = pd.read_csv(cumcounts + "/" + csv, sep=",").drop(columns=["length", "bases"])
        dfcum = dfrun.reset_index(drop=True).set_index("target")
        df_j3 = pd.DataFrame({
            "timestep": dfcum.loc["TP53_D2+J3", "timestep"].values,
            "mean": dfcum.loc["TP53_D2+J3", "mean"].values - dfcum.loc["TP53_D2_only", "mean"].values,
        })
        df_j3["target"] = "TP53_J3"
        dfcum.drop(index=["TP53_D2+J3"], inplace=True)
        dfcum.reset_index(inplace=True, drop=False)
        dfcum = pd.concat([dfcum, df_j3], ignore_index=True)
        dfcum["sample"] = csv.split("_")[0]
        cumcounts_per_sample.append(dfcum[["target", "timestep", "mean", "sample"]])  #, "experiment_id", "group_size"

df1 = pd.concat(cumcounts_per_sample, ignore_index=True)
# df1.to_csv("/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/cum_read_counts_1hr_per_target.csv",
#            index=False)

dfsup1000 = (df1[df1["mean"] >= 1000]
             .groupby(["sample", "target"])
             .agg(mean_cov=("mean", "min"), timestep=("timestep", "min"))
             .reset_index()
             .sort_values(by=["timestep"]))

dfsup1000 = (dfsup1000
             .merge(df_samples[["sample", "group_size", "experiment_id"]], left_on="sample", right_on="sample", how="left")
             .drop_duplicates())

dfsup1000.to_csv("/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/targets_above_1000x_per_sample.csv", index=False)

df1Kcount = (dfsup1000
             .groupby(["group_size", "experiment_id", "target", "timestep"])
             .count()["sample"]
             .reset_index(drop=False))
print(df1Kcount)
df1Kcumcount = df1Kcount.assign(cumcount=df1Kcount
                                .groupby(["group_size", "experiment_id", "target"])["sample"]
                                .cumsum())
print(df1Kcumcount)
df1Kcumcount.to_csv(
    "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/targets_above_1000x_over_time.csv", index=False)

m15_2 = df1Kcumcount[df1Kcumcount["experiment_id"] == "M15_2"]
#print(m15_2)
mwash4 = df1Kcumcount[df1Kcumcount["experiment_id"] == "Mwash4"]
print(mwash4)

experiments = dict([("M15_2", m15_2), ("Mwash4", mwash4)])
experiments = dict([(exp, df1Kcumcount[df1Kcumcount["experiment_id"] == exp]) for exp in df1Kcumcount["experiment_id"]])

for exp_name, exp_data in experiments.items():
    fig1, ax1 = plt.subplots(len(relabels), 1,   # as many rows as targets
                             figsize=(5, len(relabels) * 1.95),
                             sharex=True)
    for a, lab in enumerate(relabels.keys()):
        sns.lineplot(data=exp_data[exp_data["target"] == lab], x="timestep", y="cumcount", marker="o", ax=ax1[a])
        ax1[a].set_title(f'Target: {relabels[lab]}')
        ax1[a].set_xlabel('Time (minutes)')
        ax1[a].set_ylabel('')
        ax1[a].set_yticks(range(0, exp_data["cumcount"].max() + 1 + 1, 1))
        ax1[a].grid(True)
    plt.tight_layout()
    fig1.suptitle(f'Number of samples with coverage >= 1000x over time (Experiment {exp_name}: {exp_data["cumcount"].max()} samples)', y=1.02)
    ax1[len(relabels) // 2].set_ylabel('Number of samples with coverage >= 1000x')
    plt.savefig(f"/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/targets_above_1000x_over_time_{exp_name}.{fig_format}",
                bbox_inches='tight')
    plt.close()