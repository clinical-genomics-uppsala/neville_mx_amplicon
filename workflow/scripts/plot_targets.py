import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import pandas as pd

dfbed = pd.read_csv("/home/camille/ampliconthemato/amplic_ont_hemato/data/primer_data/amplicons.bed",
                    sep='\t',
                    header=None,
                    names=["chrom", "start", "end", "name", "pool"]
                    )

nchrom = len(set(dfbed["chrom"]))
fig, ax = plt.subplots(1, 1, figsize=(30, 10))
fcolors = ['b', 'g', 'r', 'c', 'm', 'y']

ax.set_xlim(left=7667400, right=7691400)
for i, target in enumerate(dfbed[dfbed["chrom"] == "chr17"].itertuples()):
    print(i, target)
    rect = ax.axvspan(xmin=target.start, xmax=target.end,
                      ymin=0.1 + 0.08 * i, ymax=0.15 + 0.08 * i,
                      facecolor=fcolors[target.pool])
    ax.text(rect.get_center()[0], rect.get_center()[1], target.name,
            fontsize=16, horizontalalignment='center', verticalalignment='center')
ax.yaxis.set_tick_params(labelleft=False)
ax.set_yticks([])

legend_elements = [Patch(facecolor='b', label="All pools"),
                   Patch(facecolor='g', label="Pool 1"),
                   Patch(facecolor='r', label="Pool 2"),
                   Patch(facecolor='c', label="Pool 3")
                   ]
ax.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(1, 0.5), fontsize=20)

plt.savefig("/home/camille/ampliconthemato/pipeline_pool_amplicon/images/TP53_mqc.png")
plt.show()
