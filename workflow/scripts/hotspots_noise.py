__author__ = "Camille Clouard"
__copyright__ = "Copyright 2026, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL3"

import pandas as pd
import json
import argparse

def parse_args():
    parser = argparse.ArgumentParser(
        description="Process hotspots noise and artifacts VAF data.",
        fromfile_prefix_chars='@'
    )
    parser.add_argument(
        "--callers",
        nargs="+",
        help="List of variant callers to process."
    )
    parser.add_argument("--hotspots_tab", help="Path to the Hotspots CSV file.")
    parser.add_argument("--background_hotpots", help="Path to the output background Hotspots CSV file.")
    parser.add_argument("--artifacts_hotpots", help="Path to the output artifacts Hotspots CSV file.")
    parser.add_argument("--background_tsv", nargs="+", help="Paths to the background TSV files.")
    parser.add_argument("--background_json", nargs="+", help="Paths to the output background JSON files.")
    parser.add_argument("--artifacts_tsv", help="Path to the artifacts TSV file.")
    parser.add_argument("--artifacts_json", help="Path to the output artifact JSON file.")
    
    return parser.parse_args()

args = parse_args()

callers = args.callers
hotspots_tab = args.hotspots_tab
background_hotpots = args.background_hotpots
artifacts_hotpots = args.artifacts_hotpots
background_tsv = args.background_tsv
background_json = args.background_json
artifacts_tsv = args.artifacts_tsv
artifacts_json = args.artifacts_json

df_hotspots = pd.read_csv(hotspots_tab, sep=",").drop(columns="Sample").set_index("Position", drop=True).drop_duplicates()
print(df_hotspots)
df_hotspots_background = df_hotspots.copy()

# Add background data for each caller to the hotspots dataframe
for caller, tsv in zip(callers, background_tsv):
    df_background = pd.read_csv(tsv, sep="\t")
    df_background.to_json(background_json[callers.index(caller)], orient="records", indent=4)
    df_background["Position"] = "chr" + df_background["Chr"].astype(str) + ":" + df_background["Pos"].astype(str)
    df_background.set_index("Position", inplace=True, drop=True)
    df_background.drop(columns=["Chr", "Pos"], inplace=True)
    df_background.rename(columns={"Median_or_Mean": f"Median or Mean VAF observed with {caller}",
                                  "SD": f"Standard deviation in VAF observed with {caller}",
                                  "NrObs": f"Number of samples with alt allele observed with {caller}"},
                         inplace=True)
    df_hotspots_background = df_hotspots_background.join(df_background, how="left", rsuffix=f"_{caller}")
print(df_hotspots_background)
df_hotspots_background.fillna(0, inplace=True)
df_hotspots_background.to_csv(background_hotpots, index=True)

# Add artifact data to the hotspots dataframe
df_hotspots_artifacts = df_hotspots.copy()
df_artifacts = pd.read_csv(artifacts_tsv, sep="\t")
df_artifacts.to_json(artifacts_json, orient="records", indent=4)
df_artifacts["Position"] = df_artifacts["Chromosome"].astype(str) + ":" + df_artifacts["pos"].astype(str)
df_artifacts.set_index("Position", inplace=True, drop=True)
df_artifacts.drop(columns=["Chromosome", "pos"], inplace=True)
df_artifacts.rename(columns={"median_MAF": f"Median VAF observed in artifacts",
                             "sd_MAF": f"Standard deviation in VAF observed in artifacts"},
                    inplace=True)
df_hotspots_artifacts = df_hotspots_artifacts.join(df_artifacts, how="left", rsuffix="_artifacts")
print(f"Artifact data from {artifacts_tsv}:", df_artifacts)
df_hotspots_artifacts.fillna(0, inplace=True)
df_hotspots_artifacts.to_csv(artifacts_hotpots, index=True)
