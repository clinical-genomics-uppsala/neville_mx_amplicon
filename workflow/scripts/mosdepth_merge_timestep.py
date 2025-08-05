import os
import pandas as pd

cols = ["length", "bases", "mean", "target"]
df = pd.concat([pd.read_csv(summary, sep='\t')
               .assign(target=os.path.basename(summary).replace(".mosdepth.summary.txt", "")) \
                for summary in list(snakemake.input)
                ])
print(df)
try:
    df = df[df["chrom"] == "total_region"].loc[:, cols]
except:  # this might not be necessary
    df = pd.DataFrame([0.0, 0.0, 0.0, 0.0], columns=cols)
print(df)
try:
    df.set_index("target", inplace=True)
except AttributeError:  # if pd.Series is returned
    df.rename("target", inplace=True)
df.to_csv(snakemake.output.csv, index=True)