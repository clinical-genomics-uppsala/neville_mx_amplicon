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
from scipy.stats import pearsonr
import os

fig_format = "svg"
date = "260622"
cohort = "CGU_2024_05"
vaf_csv = f"/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/VAF_TM_ONT_{date}_{cohort}.csv"
diff_ratio_csv = f"/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/VAF_TM_ONT_diff_ratio_{cohort}.csv"
plot1 = f"/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/VAF_diff-ratio_{cohort}_all.{fig_format}"
plot2 = f"/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/VAF_diff-ratio_{cohort}_Flongle.{fig_format}"
plot3 = f"/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/VAF-diff-ratio_{cohort}_MinION.{fig_format}"

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
print(f"Samples in the dataset: {df1['sample'].unique()} ({df1['sample'].unique().size} samples)")
df1["diff_ratio"] = abs(df1["VAF_SR_TM"] - df1["VAF_LR_ONT"]) / df1["VAF_SR_TM"]
df1["diff"] = abs(df1["VAF_SR_TM"] - df1["VAF_LR_ONT"])
df1.to_csv(diff_ratio_csv, index=False)
# df1["log_diff_ratio"] = df1["diff_ratio"].apply(lambda x: np.log10(x) if x > 0 else -np.log10(-x))
df1["log_diff_ratio"] = df1["diff_ratio"].apply(lambda x: np.log(1-x) if x <= 0 else -np.log(1-x))

df1.dropna(subset=["VAF_SR_TM", "VAF_LR_ONT"], inplace=True)
data_points_counts = df1["flowcell"].value_counts().to_dict()
# r2_vaf = r2_score(df1["VAF_SR_TM"], df1["VAF_LR_ONT"])
r2_vaf = pearsonr(df1["VAF_SR_TM"], df1["VAF_LR_ONT"]).statistic ** 2
p_val = pearsonr(df1["VAF_SR_TM"], df1["VAF_LR_ONT"]).pvalue

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

plt.annotate(text=f'R²={r2_vaf:.3f}\np={p_val:.3e}',
             xy=(100, 100), xytext=(-8, +0.90), xycoords='axes fraction',
             color='k', fontsize=16)
try:
    legend_elements = [Line2D([0], [0], marker=marker_map['Flongle'], color='w',
                              label=f"Flongle (n={data_points_counts['Flongle']})",
                              markerfacecolor='k', markersize=marker_size),
                       Line2D([0], [0], marker=marker_map['MinION'], color='w',
                              label=f"MinION (n={data_points_counts['MinION']})",
                              markerfacecolor='k', markersize=marker_size),
                       ]
except (KeyError):
    legend_elements = []
    print("No flowcell data available for this cohort, skipping the plots.")

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


# Scatterplot for Flongle data only

fig2, ax2 = plt.subplots(1, 1, figsize=(10, 9))
try:
    df_flongle = df1[df1["flowcell"] == "Flongle"]
    data_points_counts = df_flongle["flowcell"].value_counts().to_dict()
    # r2_vaf_flongle = r2_score(df_flongle["VAF_SR_TM"], df_flongle["VAF_LR_ONT"])
    r2_vaf_flongle = pearsonr(df_flongle["VAF_SR_TM"], df_flongle["VAF_LR_ONT"]).statistic ** 2
    p_val_flongle = pearsonr(df_flongle["VAF_SR_TM"], df_flongle["VAF_LR_ONT"]).pvalue
    g2 = sns.scatterplot(data=df_flongle, x="VAF_SR_TM", y="VAF_LR_ONT", hue="diff_ratio",
                         ax=ax2,
                         style="flowcell", markers=marker_map,
                         palette=cpalette, alpha=0.9, legend=False)
    ax2.plot([0, 100], [0, 100], '--k', label="1:1")
    cax2 = colorbar.make_axes(ax2, orientation='vertical', pad=0.05)
    colorbar.Colorbar(cax2[0], cmap=cpalette,
                      values=np.linspace(df_flongle["diff_ratio"].min(), df_flongle["diff_ratio"].max(), 10),
                      ).set_label(label="Magnitude of the difference in VAF", size=18)
    cax2[0].tick_params(labelsize=16)  # size of the tick labels on the colorbar
    ax2.set_xlabel("VAF% NGS", fontsize=18)
    ax2.set_ylabel("VAF% ONT", fontsize=18)
    ax2.tick_params(axis='both', which='major', labelsize=18)
    ax2.set_title("Concordance of the VAF for variants called in NGS vs. in ONT data", fontsize=20)
    plt.annotate(text=f'R²={r2_vaf_flongle:.3f}\np={p_val_flongle:.3e}',
                 xy=(100, 100), xytext=(-8, +0.90), xycoords='axes fraction',
                 color='k', fontsize=16)
    legend_elements = [Line2D([0], [0], marker=marker_map['Flongle'], color='w',
                              label=f"Flongle (n={data_points_counts['Flongle']})",
                              markerfacecolor='k', markersize=marker_size),
                       # Line2D([0], [0], marker=marker_map['MinION'], color='w',
                       #        label=f"MinION (n={data_points_counts['MinION']})",
                       #        markerfacecolor='k', markersize=marker_size),
                       ]
    ax2.legend(handles=legend_elements, title="Flowcell type", title_fontsize=14, fontsize=14, loc='upper left')
except (KeyError, ValueError):
    print("No Flongle data available for this cohort, skipping the Flongle-only plot.")
    plt.annotate(text="No Flongle data available\nfor this cohort",
                 xy=(100, 100), xytext=(-8, +0.90), xycoords='axes fraction',
                 color='k', fontsize=16)

plt.savefig(plot2, bbox_inches='tight', dpi=600)

# Scatterplot for MinION data only

fig3, ax3 = plt.subplots(1, 1, figsize=(10, 9))
df_minion = df1[df1["flowcell"] == "MinION"]
data_points_counts = df_minion["flowcell"].value_counts().to_dict()
# r2_vaf_minion = r2_score(df_minion["VAF_SR_TM"], df_minion["VAF_LR_ONT"])
r2_vaf_minion = pearsonr(df_minion["VAF_SR_TM"], df_minion["VAF_LR_ONT"]).statistic ** 2
p_val_minion = pearsonr(df_minion["VAF_SR_TM"], df_minion["VAF_LR_ONT"]).pvalue
g3 = sns.scatterplot(data=df_minion, x="VAF_SR_TM", y="VAF_LR_ONT", hue="diff_ratio",
                     ax=ax3,
                     style="flowcell", markers=marker_map,
                     palette=cpalette, alpha=0.9, legend=False)
ax3.plot([0, 100], [0, 100], '--k', label="1:1")
cax3 = colorbar.make_axes(ax3, orientation='vertical', pad=0.05)
colorbar.Colorbar(cax3[0], cmap=cpalette,
                  values=np.linspace(df_minion["diff_ratio"].min(), df_minion["diff_ratio"].max(), 10),
                  ).set_label(label="Magnitude of the difference in VAF", size=18)
cax3[0].tick_params(labelsize=16)  # size of the tick labels on the colorbar
ax3.set_xlabel("VAF% NGS", fontsize=18)
ax3.set_ylabel("VAF% ONT", fontsize=18)
ax3.tick_params(axis='both', which='major', labelsize=18)
ax3.set_title("Concordance of the VAF for variants called in NGS vs. in ONT data", fontsize=20)
plt.annotate(text=f'R²={r2_vaf_minion:.3f}\np={p_val_minion:.3e}',
             xy=(100, 100), xytext=(-8, +0.90), xycoords='axes fraction',
             color='k', fontsize=16)
legend_elements = [  # Line2D([0], [0], marker=marker_map['Flongle'], color='w',
                   #        label=f"Flongle (n={data_points_counts['Flongle']})",
                   #        markerfacecolor='k', markersize=marker_size),
                   Line2D([0], [0], marker=marker_map['MinION'], color='w',
                          label=f"MinION (n={data_points_counts['MinION']})",
                          markerfacecolor='k', markersize=marker_size),
                   ]
ax3.legend(handles=legend_elements, title="Flowcell type", title_fontsize=14, fontsize=14, loc='upper left')
plt.savefig(plot3, bbox_inches='tight', dpi=600)
