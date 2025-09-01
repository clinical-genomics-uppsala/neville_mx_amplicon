import os
import pandas as pd

"""
For each timestep, this script aggregates the average coverage for each target from mosdepth summary files and writes the values to a CSV file.
"""

os.makedirs(snakemake.output.outdir, exist_ok=True)

summaries = []
for nbamdir in os.listdir(snakemake.input.indir):
    for target in snakemake.config.get("amplicons") + snakemake.config.get("extra_regions"):
        summary = os.path.join(snakemake.input.indir, nbamdir, f"{target}.mosdepth.summary.txt")
        if os.path.exists(summary):
            summaries.append(summary)
        else:
            print(f"Warning: {summary} does not exist. Skipping.")

    cols = ["length", "bases", "mean", "target"]
    df = pd.concat([pd.read_csv(summary, sep='\t')
                   .assign(target=os.path.basename(summary).replace(".mosdepth.summary.txt", ""))
                    for summary in summaries
                    ])
    print(df)
    try:
        df = df[df["chrom"] == "total_region"].loc[:, cols]
    except KeyError:  # this might not be necessary
        df = pd.DataFrame([0.0, 0.0, 0.0, 0.0], columns=cols)
    print(df)
    try:
        df.set_index("target", inplace=True)
    except AttributeError:  # if pd.Series is returned
        df.rename("target", inplace=True)
    df.to_csv(
        os.path.join(snakemake.output.outdir, f"timestep{os.path.basename(nbamdir)}_coverage_per_amplicon.csv"),
        index=True
    )
