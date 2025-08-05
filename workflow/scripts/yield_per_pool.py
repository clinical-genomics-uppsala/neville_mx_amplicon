import pandas as pd

cols = ["pool", "target", "reads_counts"]
df = pd.read_csv(snakemake.input.csv, sep=',', index_col="target")
counts = dict([(target, round(df.loc[target, "mean"]))
               for target in df.index
               if
               target.replace(f"{snakemake.wildcards.sample}_{snakemake.wildcards.type}_", "")
               in snakemake.config.get("pools")[int(snakemake.wildcards.pooln)]
               # dirty fix to not include the sample's name in the target's label
               ])
# d_and_j = [target for target in counts.keys() if target.find("+J") >= 0]
# d_only = [(target, round(df.loc[target, "mean"]))
#           for target in df.index
#           if (target.find("TP53_D") >= 0 and target.find("_only") >= 0)
#           ]
# if d_and_j:
#     counts["TP53_JX"] = counts[d_and_j[0]] - d_only[0][1]
#     del counts[d_and_j[0]]
counts["total"] = sum(counts.values())
df_p = pd.DataFrame.from_dict(counts, orient="index", columns=["reads_counts"])
df_p.index.name = f"Pool {snakemake.wildcards.pooln}"
df_p = df_p.assign(pct_reads=round(df_p["reads_counts"] * 100 / counts["total"], ndigits=1))
df_p.to_csv(snakemake.output.csv, index=True)
