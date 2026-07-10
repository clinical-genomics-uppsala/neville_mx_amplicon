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
import re

img_extension = "svg"
timestep = 60  # in minutes
summary_csv = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/list_exp_id_sample_id.csv"
cumcounts = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/cum_read_counts"
tabcounts1 = f"/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/tab_read_counts_{timestep}minutes.csv"
boxplot1 = f"/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/boxplot_read_counts_{timestep}minutes.{img_extension}"
swarmplot1 = f"/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/swarmplot_read_counts_{timestep}minutes.{img_extension}"
boxplot4 = f"/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/boxplot_read_counts_30+60minutes.{img_extension}"

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

summary_df = pd.read_csv(summary_csv, sep=",")
print(summary_df)

# Read cumulative counts after 1 hour

cumcounts_Xminutes = []
for csv in os.listdir(cumcounts):
    if csv[-4:] == ".csv":
        print(csv)
        try:
            sample = re.compile(r'D\d{2}\-\d{5}(\-\d)*').search(csv).group()
            experiment = csv.split(f"_{sample}")[0]
            flowcell = csv.split(sample)[1].split("_")[1]
            print(sample, experiment, flowcell)
        except AttributeError:
            print(f"Could not parse sample and flowcell from {csv}, skipping.")
            continue
        if flowcell == "MinION" and experiment in summary_df["experiment_id"].values:
            dfrun = pd.read_csv(cumcounts + "/" + csv, sep=",")
            # Wash5 and M22_run2 have 60-minute timesteps
            if dfrun["timestep"].min() == 60 and experiment not in ["Wash5", "M22_run2"]:
                dfrun["timestep"] = dfrun["timestep"] // 6
            df1hour = dfrun[dfrun["timestep"] == timestep].reset_index(drop=True).set_index("target")
            print(dfrun)
            print(df1hour)
            if not df1hour.empty:
                counts_j3 = {"target": "TP53_J3", "mean": df1hour.loc["TP53_D2+J3", "mean"] - df1hour.loc["TP53_D2_only", "mean"]}
                df1hour.drop(index=["TP53_D2+J3"], inplace=True)
                df1hour.reset_index(inplace=True, drop=False)
                df1hour.loc[len(df1hour)] = counts_j3
                df1hour["sample"] = sample
                df1hour["experiment"] = experiment
                # fetch group_size from summary csv
                cumcounts_Xminutes.append(df1hour[["target", "mean", "sample", "experiment"]])

print(f"Read cumulative counts after 1 hour sequencing on MinION flowcells for {len(cumcounts_Xminutes)} samples.")
print(cumcounts_Xminutes)
df1h = pd.concat(cumcounts_Xminutes, ignore_index=True)
exp_to_group = dict(
    (t.experiment_id, t.group_size)
    for t in summary_df[["experiment_id", "group_size"]].itertuples()
)
df1h["group_size"] = df1h["experiment"].replace(exp_to_group)
df1h["target"] = df1h["target"].replace(relabels)
print(df1h)
df1h.to_csv(tabcounts1, index=False)
# print("\nTargets: ", set(df1h["target"]))

sns.set_style("whitegrid")
resampled_palette = mpl.colormaps['plasma'].resampled(8)
colors = resampled_palette(np.linspace(0, 0.7, df1h["group_size"].nunique()))
# exclude interval [0.6, 1) to avoid picking yellow shades

fig1h, ax1h = plt.subplots(1, 1, figsize=(12, 8))
leglabels = []
leghandles = []
for t in df1h["group_size"]:
    try:
        t = int(t)
    except ValueError:
        print(f"Could not convert {t} to int, skipping.")
        continue
for i, group in enumerate(sorted(df1h["group_size"].astype(int).unique())):
    sns.stripplot(x="target",
                  y="mean",
                  size=8,
                  data=df1h[df1h["group_size"] == group],
                  color=colors[i],
                  alpha=0.8,
                  ax=ax1h,
                  legend=False
                  )
    leglabels.append(f"{group} sample(s)")
    leghandles.append(Line2D([], [], color="white",
                             marker='o', markerfacecolor=colors[i], markersize=10,))
ax1h.set_xticklabels(ax1h.get_xticklabels(),
                     rotation=30,
                     ha='right',
                     )
# labs1h = ax1h.get_xticklabels()
# newlabs1h = [relabels[lab.get_text()] for lab in labs1h]
# ax1h.set_xticklabels(newlabs1h)
ax1h.set_ylim(bottom=0, top=df1h["mean"].max() + 100)
plt.legend(handles=leghandles, labels=leglabels, title="Group size",
           title_fontsize=14, fontsize=12, loc='upper right')
for target in set(df1h["target"]):
    plt.axvline(x=list(set(df1h["target"])).index(target) - 0.0,
                linestyle='--', color='lightgrey', alpha=0.5)
plt.title(f"Estimated read counts after {timestep} minutes sequencing on MinION flowcells (n={len(df1h['sample'])})",
          fontsize=18)
ax1h.set_xlabel("Amplicon", fontsize=16)
ax1h.set_ylabel("Read counts", fontsize=16)
plt.savefig(swarmplot1, bbox_inches='tight')
plt.show()
sns.set_style("white")


_n_rows = len(df1h["group_size"].unique())
boxfig1h, boxax1h = plt.subplots(_n_rows, 1, figsize=(12, 8 * _n_rows))  # , sharex=True)
leglabels = []
leghandles = []
for i, group in enumerate(sorted(df1h["group_size"].unique())):
    sns.boxplot(x="target",
                y="mean",
                data=df1h[df1h["group_size"] == group],
                # color=colors[i],
                ax=boxax1h[i],
                legend=False
                )
    leglabels.append(f"{group} sample(s)")
    leghandles.append(Line2D([], [], color="white",
                             marker='o', markerfacecolor=colors[i], markersize=10,))
    boxax1h[i].set_xticklabels(
        df1h[df1h["group_size"] == group]["target"].unique(),  # boxax1h[i].get_xticklabels(),
        rotation=30,
        ha='right',
    )
    print(df1h[df1h["group_size"] == group])
    boxax1h[i].set_ylim(bottom=0, top=df1h[df1h["group_size"] == group]["mean"].max() + 100)
    boxax1h[i].set_xlabel("Amplicon", fontsize=16)
    boxax1h[i].set_ylabel("Read counts", fontsize=16)
    boxax1h[i].axhline(y=1000, linestyle='--', color='r', label='1000')
    # boxax1h[i].legend(handles=[f"{group} sample(s)"], labels=[Line2D([], [], color="white",
    #                                                                  marker='o', markerfacecolor=colors[i], markersize=10,)],
    #                   title="Group size",
    #                   title_fontsize=14, fontsize=12, loc='upper right')
    for target in set(df1h["target"]):
        boxax1h[i].axvline(
            x=list(set(df1h["target"])).index(target) - 0.0,
            linestyle='--', color='lightgrey', alpha=0.5
        )
    boxax1h[i].set_title(f"(group size={group}, n={len(df1h[df1h['group_size'] == group]['sample']) // len(relabels)} samples)",
                         fontsize=18)
plt.suptitle(f"Estimated read counts after {timestep} minutes sequencing on MinION flowcell", fontsize=20, y=0.92)
# plt.tight_layout()
plt.savefig(boxplot1, bbox_inches='tight')
plt.show()
sns.set_style("white")

tsteps = [60, 120]
try:
    df2 = pd.read_csv(f"/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/tab_read_counts_60minutes.csv", sep=",")
    df3 = pd.read_csv(f"/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/tab_read_counts_120minutes.csv", sep=",")
    df2["timestep"] = 60
    df3["timestep"] = 120
    df4 = pd.concat([df2, df3], ignore_index=True)
    df4["target"] = df4["target"].replace(relabels)
    print(df4)
    min_reads = df4["mean"].min()
    print(df4[df4["mean"] <= min_reads + 200])
    print(df4.groupby(["target", "timestep"]).describe())
    (df4.groupby(["target", "timestep"])
     .describe()
     .to_csv("/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/4+24hours_describe.csv", index=False))

    sns.set_theme(style="ticks")
    custom_palette = {60: "darkkhaki", 120: "royalblue"}
    fig4, ax4 = plt.subplots(_n_rows, 1, figsize=(12, 8 * _n_rows))  # , sharex=True)

    for i, group in enumerate(sorted(df1h["group_size"].unique())):
        sns.boxplot(x="target",
                    y="mean",
                    hue="timestep",
                    data=df4[df4["group_size"] == group],
                    palette=custom_palette,
                    width=0.8,
                    linewidth=2.0,
                    showmeans=True,
                    showfliers=False,
                    # fliersize=1,
                    ax=ax4[i],
                    )
        hds, labs = ax4[i].get_legend_handles_labels()
        hours = []
        for hd, lab in zip(hds, labs):
            lab = f"{int(lab) // 60} hours"
            hours.append(lab)
        ax4[i].legend(handles=hds, labels=hours,
                      title_fontsize=18, fontsize=18)
        ax4[i].set_xticklabels(ax4[i].get_xticklabels(),
                               rotation=30,
                               ha='right',
                               )
        ax4[i].tick_params(axis='both', which='major', labelsize=14)
        ax4[i].set_xlabel("Amplicon", fontsize=18)
        ax4[i].set_ylabel("Read counts", fontsize=18)
        # ax4[i].set_ylim(bottom=0, top=df1h[df1h["group_size"] == group]["mean"].max() + 100)
        ax4[i].axhline(y=1000, linestyle='--', color='r', label='1000')
        ax4[i].set_title(f"(group size={group}, n={len(df1h[df1h['group_size'] == group]['sample']) // len(relabels)} samples)",
                         fontsize=18)
        for target in set(df1h["target"]):
            ax4[i].axvline(x=list(set(df1h["target"])).index(target) - 0.0,
                           linestyle='--', color='lightgrey', alpha=0.5)
        # red dashed line at y=1000 reads
        # plt.annotate(text='1000 reads', xy=(0, 1000), xytext=(-0.11, 0.02), xycoords='figure fraction', color='r')
    plt.suptitle(f"Estimated read counts after 60 and 120 minutes sequencing on MinION flowcell", fontsize=20, y=0.92)
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.5)
    plt.savefig(boxplot4, bbox_inches='tight')

except FileNotFoundError:
    pass

if False:
    # Read cumulative counts after 4 hours

    cumcounts_4hours = []
    for csv in os.listdir(cumcounts):
        if csv[-4:] == ".csv":
            df24h = pd.read_csv(cumcounts + "/" + csv, sep=",")
            df4hours = df24h[df24h["timestep"] == 240].reset_index(drop=True).set_index("target")
            counts_j3 = {"target": "TP53_J3", "mean": df4hours.loc["TP53_D2+J3", "mean"] - df4hours.loc["TP53_D2_only", "mean"]}
            df4hours.drop(index=["TP53_D2+J3"], inplace=True)
            df4hours.reset_index(inplace=True, drop=False)
            df4hours.loc[len(df4hours)] = counts_j3
            df4hours["sample"] = csv.split("_")[0]
            cumcounts_4hours.append(df4hours[["target", "mean", "sample"]])

    df2 = pd.concat(cumcounts_4hours, ignore_index=True)
    print("\nTargets: ", set(df2["target"]))
    print(df2.groupby(["target"]).describe())

    plt.style.use('classic')
    fig2, ax2 = plt.subplots(1, 1, figsize=(14, 6))
    df2.groupby(["target"]).boxplot(rot=90, fontsize=12, ax=ax2, subplots=False)
    labs2 = ax2.get_xticklabels()
    newlabs2 = []
    for lab in labs2:
        newlabs2.append(relabels[lab.get_text().strip('()mean').replace(",", "").strip()])
    ax2.set_xticklabels(newlabs2)
    ax2.set_ylabel("Read counts")
    ax2.set_title("Estimated read counts after 4 hours sequencing")
    plt.savefig(plot2, bbox_inches='tight')

    # Read cumulative counts after 24 hours
    cumcounts_24hours = []
    for csv in os.listdir(cumcounts):
        if csv[-4:] == ".csv":
            df24h = pd.read_csv(cumcounts + "/" + csv, sep=",")
            df24hours = df24h[df24h["timestep"] == 1440].reset_index(drop=True).set_index("target")
            counts_j3 = {"target": "TP53_J3", "mean": df24hours.loc["TP53_D2+J3", "mean"] - df24hours.loc["TP53_D2_only", "mean"]}
            df24hours.drop(index=["TP53_D2+J3"], inplace=True)
            df24hours.reset_index(inplace=True, drop=False)
            df24hours.loc[len(df24hours)] = counts_j3
            df24hours["sample"] = csv.split("_")[0]
            cumcounts_24hours.append(df24hours[["target", "mean", "sample"]])

    df3 = pd.concat(cumcounts_24hours, ignore_index=True)
    print(df3.groupby(["target"]).describe())
    min_24h = df3["mean"].min()
    print("Min read counts after 24h: ", df3[df3["mean"] <= min_24h + 1000])
    print(df3.sort_values(by="mean", ascending=True))
    print(df3.sort_values(by="mean", ascending=False))

    plt.style.use('classic')
    fig3, ax3 = plt.subplots(1, 1, figsize=(12, 6))
    df3.groupby(["target"]).boxplot(rot=90, fontsize=12, ax=ax3, subplots=False)
    labs3 = ax2.get_xticklabels()
    newlabs3 = []
    for lab in labs3:
        newlabs3.append(lab.get_text().strip('()mean').replace(",", "").strip().replace(",", ""))
    ax3.set_xticklabels(newlabs3)
    ax3.set_ylabel("Read counts")
    ax3.set_title("Estimated read counts after 24 hours sequencing")
    plt.savefig(plot3, bbox_inches='tight')

    # Combine 4 hours and 24 hours data for boxplot and swarmplot

    df2["timestep"] = 240
    df3["timestep"] = 1440
    df4 = pd.concat([df2, df3], ignore_index=True)
    df4["target"] = df4["target"].replace(relabels)
    print(df4)
    min_reads = df4["mean"].min()
    print(df4[df4["mean"] <= min_reads + 200])
    print(df4.groupby(["target", "timestep"]).describe())
    (df4.groupby(["target", "timestep"])
     .describe()
     .to_csv("/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/4+24hours_describe.csv", index=False))

    sns.set_theme(style="ticks")
    custom_palette = {240: "darkkhaki", 1440: "royalblue"}
    fig4, ax4 = plt.subplots(1, 1, figsize=(16, 12))
    # df4.groupby(["target", "timestep"]).boxplot(rot=90, fontsize=12, ax=ax3, subplots=False)
    sns.boxplot(x="target",
                y="mean",
                hue="timestep",
                data=df4,
                palette=custom_palette,
                width=0.8,
                linewidth=2.0,
                showmeans=True,
                fliersize=1,
                ax=ax4,
                )
    hds, labs = ax4.get_legend_handles_labels()
    hours = []
    for hd, lab in zip(hds, labs):
        lab = f"{int(lab) // 60} hours"
        hours.append(lab)
    plt.legend(handles=hds, labels=hours,
               title_fontsize=18, fontsize=18)
    ax4.set_xticklabels(ax4.get_xticklabels(),
                        rotation=30,
                        ha='right',
                        )
    ax4.tick_params(axis='both', which='major', labelsize=14)
    ax4.set_xlabel("Amplicon", fontsize=18)
    ax4.set_ylabel("Read counts", fontsize=18)
    ax4.set_title("Estimated read counts after 4 and 24 hours sequencing", fontsize=18)
    ax4.set_ylim(bottom=0, top=55000)
    plt.axhline(y=1000, linestyle='--', color='r', label='1000')
    plt.annotate(text='1000 reads', xy=(0, 1000), xytext=(-0.11, 0.02), xycoords='axes fraction', color='r')
    plt.savefig(boxplot4, bbox_inches='tight')

    ax4.clear()
    sns.swarmplot(x="target",
                  y="mean",
                  hue="timestep",
                  data=df4,
                  palette=custom_palette,  # Use the custom palette
                  ax=ax4,
                  )
    ax4.set_xticklabels(ax4.get_xticklabels(),
                        rotation=30,
                        ha='right',
                        )
    # labs4 = ax2.get_xticklabels()
    # newlabs4 = []
    # for lab in labs4:
    #     newlabs4.append(relabels[lab.get_text().strip('()mean').replace(",", "").strip()])
    ax4.set_ylim(bottom=0, top=None)
    plt.axhline(y=1000, linestyle='--', color='r', label='1000 threshold')
    plt.savefig(swarmplot4, bbox_inches='tight')
