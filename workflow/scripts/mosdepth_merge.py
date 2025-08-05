import os
import pandas as pd

cols = ["length", "bases", "mean", "target"]
df = pd.concat([pd.read_csv(summary, sep='\t')
                .assign(target=os.path.basename(summary).replace(".mosdepth.summary.txt", ""))
                for summary in list(snakemake.input)
                ]).loc[3, cols]  # lines for "total_region" have index=3 in the dataframe
df.set_index("target", inplace=True)
df.to_csv(snakemake.output.csv, index=True)
