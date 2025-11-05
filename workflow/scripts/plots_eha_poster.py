import matplotlib.pyplot as plt
from matplotlib import colorbar
import matplotlib as mpl
import seaborn as sns
import pandas as pd
import numpy as np
from sklearn.metrics import r2_score
import os

vaf_csv = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/VAF_TM_ONT_251029.csv"
plot1 = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/VAF_diff-ratio.png"

multiqc_data = "/home/camille/Documents/CGU_2024_05-IDH-TP53-NPM1-nanopore/multiqc_data"

len_bins = [0, 100, 500, 1000, 2000, 2500, 3000, 3500, 4000, 5000, 6000, 7000, 8000, 500000,]

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
    'TP53_D2_only': 'TP53_D'
}

len_data = []

for seqrun in os.listdir(multiqc_data):
    if seqrun.endswith("multiqc_amplicons_data"):
        print(f"\nFetching data in {seqrun}")
        sample_name = seqrun.split('_')[0]
        dfrun = pd.read_csv(os.path.join(multiqc_data, seqrun, "pycoqc_read_len_plot_Passing_Reads.txt"), sep="\t")
        pycoqc_data = [tup.strip("()").replace(" ", "").split(",") for tup in dfrun.iloc[0].values[1:]]
        # print(df1.iloc[0].values[1:][0].strip("()").replace(" ", "").split(","))
        len_counts = pd.DataFrame(data=pycoqc_data, columns=["len", "count"], dtype=float)
        len_counts["sample"] = sample_name
        len_data.append(len_counts)

df1 = pd.concat(len_data, ignore_index=True)
intervals = pd.cut(df1["len"], bins=len_bins)
dfhist = (df1.join(intervals, rsuffix="_bin")
          .groupby(["len_bin", "sample"])["count"]
          .sum()
          .to_frame()
          .reset_index(drop=False))

fig1, ax1 = plt.subplots(1, 1, figsize=(14, 9))
sns.barplot(x="len_bin", y="count", data=dfhist,
            hue="sample", orient="v",
            ax=ax1
            )
# ax1.set_xticklabels(ax1.get_xticklabels(), rotation=45)
fmt_xlab = []
for lab in ax1.get_xticklabels():
    bounds_str = lab.get_text().strip("(]").replace(" ", "").split(",")
    fmt_xlab.append(f"{bounds_str[0]}\n"
                    f"to\n"
                    f"{bounds_str[1]}")
plt.axvline(x=3.5, linestyle='--', color='k')
plt.axvline(x=7.5, linestyle='--', color='k')
ax1.set_xticklabels(fmt_xlab)
ax1.set_xlabel("Read length (bp)", fontsize=18)
ax1.set_ylabel("Reads count", fontsize=18)
ax1.tick_params(axis='both', which='major', labelsize=16)
ax1.set_title("Length distribution of passed reads in the sequencing experiments", fontsize=20)

hds, labs = ax1.get_legend_handles_labels()
i = 1
nsamples = len(set(dfhist["sample"]))
print(nsamples)
plt.legend(title="Sample",
           handles=hds, labels=[i + 1 for i in range(nsamples)],
           title_fontsize=18, fontsize=18
           )

plt.tight_layout()
plt.show()

