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

vaf_csv = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/VAF_TM_ONT_251029.csv"
plot1 = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/VAF_diff-ratio.png"

swarmplot_minion = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/swarmplot_read_counts_1hour_minion.png"
cumcounts = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/cum_read_counts"
plot2 = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/read_counts_4hours.png"
plot3 = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/read_counts_24hours.png"
boxplot4 = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/boxplot_read_counts_4+24hours.png"
swarmplot4 = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/swarmplot_read_counts_4+24hours.png"
violinplot4 = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/violinplot_read_counts_4+24hours.png"

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

df1 = pd.read_csv(vaf_csv, sep=",")
df1["diff_ratio"] = abs(df1["VAF_SR_TM"] - df1["VAF_LR_ONT"]) / df1["VAF_SR_TM"]
df1["diff"] = abs(df1["VAF_SR_TM"] - df1["VAF_LR_ONT"])
df1.to_csv("/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/VAF_TM_ONT_diff_ratio.csv", index=False)
# df1["log_diff_ratio"] = df1["diff_ratio"].apply(lambda x: np.log10(x) if x > 0 else -np.log10(-x))
df1["log_diff_ratio"] = df1["diff_ratio"].apply(lambda x: np.log(1-x) if x <= 0 else -np.log(1-x))

df1.dropna(subset=["VAF_SR_TM", "VAF_LR_ONT"], inplace=True)
data_points_counts = df1["flowcell"].value_counts().to_dict()
r2_vaf = r2_score(df1["VAF_SR_TM"], df1["VAF_LR_ONT"])

## Scatterplot VAF NGS vs VAF ONT colored by diff_ratio

cpalette = "plasma"  # "cool"
marker_map = {'Flongle': "$\circ$",
              'MinION': 'X'}
marker_size = 14
mpl.rcParams['lines.markersize'] = marker_size
print("Font family used in plots: ", mpl.rcParams['font.family'])
font = findfont(FontProperties(family=[mpl.rcParams['font.family'][0]]))
print("Font used in plots: ", font.split('/')[-1].replace('.ttf', ''))
zoomin = False

fig1, ax1 = plt.subplots(1, 1, figsize=(10, 9))
g1 = sns.scatterplot(data=df1, x="VAF_SR_TM", y="VAF_LR_ONT", hue="diff_ratio",
                     ax=ax1,
                     style="flowcell", markers=marker_map,
                     palette=cpalette, alpha=0.9, legend=False)
ax1.plot([0, 100], [0, 100], '--k', label="1:1")
cax1 = colorbar.make_axes(ax1, orientation='vertical', pad=0.05)
colorbar.Colorbar(cax1[0], cmap=cpalette,
                  values=np.linspace(df1["diff_ratio"].min(), df1["diff_ratio"].max(), 10),
                  ).set_label(label="Magnitude of the difference in VAF", size=18)
cax1[0].tick_params(labelsize=16)  # size of the tick labels on the colorbar
ax1.set_xlabel("VAF% NGS", fontsize=18)
ax1.set_ylabel("VAF% ONT", fontsize=18)
ax1.tick_params(axis='both', which='major', labelsize=18)
ax1.set_title("Concordance of the VAF for variants called in NGS vs. in ONT data", fontsize=20)

plt.annotate(text=f'R²={r2_vaf:.3f}',
             xy=(100, 100), xytext=(-8, +0.90), xycoords='axes fraction',
             color='k', fontsize=16)

legend_elements = [Line2D([0], [0], marker="$\circ$", color='w',
                          label=f"Flongle (n={data_points_counts['Flongle']})",
                          markerfacecolor='k', markersize=marker_size),
                   Line2D([0], [0], marker='X', color='w',
                          label=f"MinION (n={data_points_counts['MinION']})",
                          markerfacecolor='k', markersize=marker_size),
                   ]
ax1.legend(handles=legend_elements, title="Flowcell type", title_fontsize=14, fontsize=14, loc='upper left')

if zoomin:
    # inset Axes....
    x1, x2, y1, y2 = -0.5, 22, -0.5, 22  # subregion of the original image
    axins = ax1.inset_axes(
        [0.52, 0.06, 0.42, 0.42],
        xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])
    g1zoom = sns.scatterplot(data=df1, x="VAF_SR_TM", y="VAF_LR_ONT", hue="diff_ratio",
                             ax=axins,
                             style="flowcell", markers=marker_map,
                             palette=cpalette, alpha=0.9, legend=False)
    axins.plot([0, 22], [0, 22], '--k')
    axins.set_title("Zoomed in on VAF 0-20%", fontsize=12)
    ax1.indicate_inset_zoom(axins, edgecolor="black")

plt.savefig(plot1, bbox_inches='tight', dpi=600)


## Read cumulative counts after 1 hour

cumcounts_1hour = []
for csv in os.listdir(cumcounts):
    if csv[-4:] == ".csv" and csv.split("_")[1] == "MinION" :
        dfrun = pd.read_csv(cumcounts + "/" + csv, sep=",")
        df1hour = dfrun[dfrun["timestep"] == 60].reset_index(drop=True).set_index("target")
        counts_j3 = {"target": "TP53_J3", "mean": df1hour.loc["TP53_D2+J3", "mean"] - df1hour.loc["TP53_D2_only", "mean"]}
        df1hour.drop(index=["TP53_D2+J3"], inplace=True)
        df1hour.reset_index(inplace=True, drop=False)
        df1hour.loc[len(df1hour)] = counts_j3
        df1hour["sample"] = csv.split("_")[0]
        cumcounts_1hour.append(df1hour[["target", "mean", "sample"]])

df1h = pd.concat(cumcounts_1hour, ignore_index=True)
df1h = df1h.merge(df1[df1["flowcell"] == "MinION"][["sample", "group_size"]].drop_duplicates(),
                  how='right',
                  on="sample")
print(df1h)
print("\nTargets: ", set(df1h["target"]))
print(df1h.groupby(["target"]).describe())

sns.set_style("whitegrid")
resampled_palette = mpl.colormaps['plasma'].resampled(8)
colors = resampled_palette(np.linspace(0, 0.7, df1h["group_size"].nunique()))
# exclude interval [0.6, 1) to avoid picking yellow shades
print(colors)
fig1h, ax1h = plt.subplots(1, 1, figsize=(12, 8))
leglabels = []
leghandles = []
for i, group in enumerate(sorted(df1h["group_size"].unique())):
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
labs1h = ax1h.get_xticklabels()
newlabs1h = [relabels[lab.get_text()] for lab in labs1h]
ax1h.set_xticklabels(newlabs1h)
ax1h.set_ylim(bottom=0, top=df1h["mean"].max() + 100)
plt.legend(handles=leghandles, labels=leglabels, title="Group size",
           title_fontsize=14, fontsize=12, loc='upper right')
for target in set(df1h["target"]):
    plt.axvline(x=list(set(df1h["target"])).index(target) - 0.0,
                linestyle='--', color='lightgrey', alpha=0.5)
plt.title(f"Estimated read counts after 1 hour sequencing on MinION flowcells (n={len(df1h['sample'].unique())})",
          fontsize=18)
ax1h.set_xlabel("Amplicon", fontsize=16)
ax1h.set_ylabel("Read counts", fontsize=16)
plt.savefig(swarmplot_minion, bbox_inches='tight')
plt.show()
sns.set_style("white")

if False:
    ## Read cumulative counts after 4 hours

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


    ## Read cumulative counts after 24 hours
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


    ## Combine 4 hours and 24 hours data for boxplot and swarmplot

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
