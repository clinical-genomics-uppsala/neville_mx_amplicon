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

# Scatterplot VAF NGS vs VAF ONT colored by diff_ratio

cpalette = "plasma"  # "cool"
marker_map = {'Flongle': r"$\circ$",
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

legend_elements = [Line2D([0], [0], marker=marker_map['Flongle'], color='w',
                          label=f"Flongle (n={data_points_counts['Flongle']})",
                          markerfacecolor='k', markersize=marker_size),
                   Line2D([0], [0], marker=marker_map['MinION'], color='w',
                          label=f"MinION (n={data_points_counts['MinION']})",
                          markerfacecolor='k', markersize=marker_size),
                   ]
ax1.legend(handles=legend_elements, title="Flowcell type", title_fontsize=14, fontsize=14, loc='upper left')

if zoomin:
    # inset Axes
    x1, x2, y1, y2 = -0.5, 22, -0.5, 22  # subregion of the original image
    axins = ax1.inset_axes(
        (0.52, 0.06, 0.42, 0.42),
        xlim=(x1, x2), ylim=(y1, y2), xticklabels=[], yticklabels=[])
    g1zoom = sns.scatterplot(data=df1, x="VAF_SR_TM", y="VAF_LR_ONT", hue="diff_ratio",
                             ax=axins,
                             style="flowcell", markers=marker_map,
                             palette=cpalette, alpha=0.9, legend=False)
    axins.plot([0, 22], [0, 22], '--k')
    axins.set_title("Zoomed in on VAF 0-20%", fontsize=12)
    ax1.indicate_inset_zoom(axins, edgecolor="black")

plt.savefig(plot1, bbox_inches='tight', dpi=600)
