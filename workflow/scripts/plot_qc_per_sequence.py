"""
The script assumes that MultiQC data is available and that sequali reports are part of those data.
It generates per-sequence QC plots for each sample based on the MultiQC data.
"""

import matplotlib.pyplot as plt
from matplotlib import colorbar
from matplotlib.lines import Line2D
from matplotlib.font_manager import findfont, FontProperties
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.metrics import r2_score
import os

mqc_data = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/multiqc_data"
qscores = "sequali_per_sequence_quality_scores_plot.txt"

# f"{mqc_data}/D24-00846_MinION_T_multiqc_amplicons_data/"
qsequence_data = {}

for sample_flowcell in os.listdir(mqc_data):
    sample = sample_flowcell.split("_")[0]
    flowcell = sample_flowcell.split("_")[1]
    seq_file = f"{mqc_data}/{sample}_{flowcell}_T_multiqc_amplicons_data/{qscores}"

    if os.path.exists(seq_file):
        df = pd.read_csv(seq_file, sep="\t").iloc[0][1:].str.split(pat=", ", expand=True)
        df[0] = df[0].str.strip('(').astype(int)
        df[1] = df[1].str.strip(')').astype(int)
        df.columns = ['Q-Score', 'Counts']
        df.set_index('Q-Score', inplace=True)
        qsequence_data[(flowcell, sample)] = df

df_qscores = pd.concat(qsequence_data.values(),
                       keys=qsequence_data.keys(),
                       names=['Flowcell', 'Sample', 'Q-Score']).reset_index()
tot_reads = df_qscores.groupby(['Flowcell', 'Sample']).agg(total_reads=("Counts", "sum"))
print(tot_reads.to_dict())
print(df_qscores)
df_pct_qscores = df_qscores.join(tot_reads, on=['Flowcell', 'Sample'])
df_pct_qscores['Pct_Counts'] = df_pct_qscores['Counts'] / df_pct_qscores['total_reads'] * 100
print(df_pct_qscores)

plt.figure(figsize=(10, 6))
sns.lineplot(data=df_pct_qscores, x='Q-Score', y='Pct_Counts',
             hue='Sample', style='Flowcell',
             alpha=0.5)
plt.title(f'Per-Sequence Quality Scores')
plt.xlabel('Mean Q-Score')
plt.ylabel('Percentage of Reads (%)')
plt.xlim(0, 50)
plt.grid(True)
plt.savefig(f"{mqc_data}/multiqc_amplicons_data_{qscores}".replace(".txt", ".png"),
            bbox_inches='tight')
plt.close()

df_agg_qscores = (df_pct_qscores.groupby(['Q-Score', "Flowcell"])
                  .agg(min_pct_counts=("Pct_Counts", "min"),
                       max_pct_counts=("Pct_Counts", "max"),
                       median_pct_counts=("Pct_Counts", "median"),
                       mean_pct_counts=("Pct_Counts", "mean"))
                  .reset_index())
df_cum_qscores = (df_agg_qscores.set_index(["Flowcell", 'Q-Score'])
                  .groupby(["Flowcell"])["mean_pct_counts"]
                  .cumsum()).reset_index()
df_q90 = df_cum_qscores[df_cum_qscores['mean_pct_counts'] <= 90].groupby('Flowcell').agg(Q90_QScore=('Q-Score', 'max'))
df_q90_qscores = (df_cum_qscores[df_cum_qscores['mean_pct_counts'] >= 90].set_index(["Flowcell", 'Q-Score'])  # [df_agg_qscores['mean_pct_counts'] <= 90]
                  .groupby(["Flowcell", 'Q-Score'])
                  .max()).reset_index()

print(df_agg_qscores)
print(df_cum_qscores)
print(df_q90_qscores)
print("Q90 Q-Scores:")
print(df_q90)

#plt.figure(figsize=(10, 6))
fig, g_agg = plt.subplots(1, sharex=True, figsize=(10, 6))
sns.lineplot(data=df_agg_qscores, x='Q-Score', y='median_pct_counts', ax=g_agg,
             hue='Flowcell',
             # style="Flowcell", markers=True, dashes=False,
             alpha=1)
sns.lineplot(data=df_agg_qscores, x='Q-Score', y='mean_pct_counts', ax=g_agg,
             hue='Flowcell', linestyle='--',
             alpha=1)
g_agg.fill_between(df_agg_qscores['Q-Score'][df_agg_qscores['Flowcell'] == 'Flongle'],
                      df_agg_qscores['min_pct_counts'][df_agg_qscores['Flowcell'] == 'Flongle'],
                      df_agg_qscores['max_pct_counts'][df_agg_qscores['Flowcell'] == 'Flongle'],
                      alpha=0.2, color=sns.color_palette()[0],
                      label='Flongle min-max range')
g_agg.fill_between(df_agg_qscores['Q-Score'][df_agg_qscores['Flowcell'] == 'MinION'],
                      df_agg_qscores['min_pct_counts'][df_agg_qscores['Flowcell'] == 'MinION'],
                      df_agg_qscores['max_pct_counts'][df_agg_qscores['Flowcell'] == 'MinION'],
                      alpha=0.2, color=sns.color_palette()[1],
                      label='MinION min-max range')
# sns.lineplot(df_cum_qscores, x='Q-Score', y='mean_pct_counts',
#             hue='Flowcell',
#             ax=g_agg[1],
#             alpha=0.5)
handles, labels = g_agg.get_legend_handles_labels()
print(handles, labels)
plt.legend(handles=handles,
           labels=['Flongle median',
                   'MinION median',
                   f'Flongle mean (Q90={df_q90.loc["Flongle", "Q90_QScore"]})',
                   f'MinION mean (Q90={df_q90.loc["MinION", "Q90_QScore"]})',
                   'Flongle min-max range',
                   'MinION min-max range'],
           title=None,
           loc='upper right')
g_agg.set_title(f'Aggregated Per-Sequence Quality Scores')
plt.xlabel('Mean Q-Score')
plt.ylabel('Percentage of Reads (%)')
plt.xlim(0, 50)
plt.grid(True)
plt.savefig(f"{mqc_data}/multiqc_amplicons_data_{qscores}".replace(".txt", "_aggregated.png"),
            bbox_inches='tight')
plt.close()