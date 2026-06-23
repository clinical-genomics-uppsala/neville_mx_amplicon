__author__ = "Camille Clouard"
__copyright__ = "Copyright 2026, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL3"

"""
This script estimates the sequencing noise in reads in a BAM file for a selection of known hotspot variants.
That is, given the genomic position of a hotspot mutation, the script will count how many reads carry the reference allele
or the alternative allele in that position,  and compute the variant allele frequency (VAF) for the alternative allele.

Columns in the input CSV file with the genomic positions of hotspot mutations should include:
- Sample: The sample identifier (e.g., Sample1, Sample2, etc.),
- Gene: The gene where the genomic position is located (e.g., TP53, NPM1, etc.),
- Position: The genomic position of the hotspot (e.g., chr1:12345),
- Transcript: The transcript associated with the gene (e.g., ENST00000269305),
- Protein_change: The protein change associated with the hotspot (e.g., p.R175H),
- Ref_allele: The reference allele at the genomic position (e.g., G).

Usage:
```bash
python workflow/scripts/hotspots_analysis.py @workflow/scripts/hotspots_analysis_args.txt
```
"""

import os
import pysam
import json
import argparse
import pandas as pd
from collections import Counter
from typing import Dict, Any, List


def load_hotspots_metadata(hspt_csv_path: str) -> Dict[str, Any]:
    """
    Load the metadata about the hotspots from the provided CSV file and return it as a dictionary.
    The CSV file is expected to have the following columns: "Position", "Gene", "Transcript", "Protein_change", "Ref_allele".
    The returned dictionary will have the genomic position as keys and a dictionary of metadata as values.

    Parameters
    ----------
    hspt_csv_path : str
        Path to the CSV file containing the hotspots metadata.

    Returns
    -------
    Dict[str, Any]
        A dictionary where keys are genomic positions (e.g. "chr13:32310000") and values are dictionaries containing the metadata for that hotspot, including:
        - "gene": gene name
        - "transcripts": comma-separated list of transcripts affected by the hotspot mutation
        - "protein_changes": comma-separated list of protein changes associated with the hotspot mutation
        - "ref_allele": reference allele at the hotspot position
    """
    df = pd.read_csv(hspt_csv_path, sep=",")
    required = {"Position", "Gene", "Transcript", "Protein_change", "Ref_allele"}
    if not required.issubset(df.columns):
        missing = required - set(df.columns)
        raise ValueError(f"Missing required columns: {missing}")

    agg = (
        df.groupby(["Position", "Gene"])
          .agg({
              "Transcript": lambda s: sorted({str(x) for x in s.dropna()}),
              "Protein_change": lambda s: sorted({str(x) for x in s.dropna()}),
              "Ref_allele": lambda s: sorted({str(x) for x in s.dropna()})
          })
          .reset_index()
    )

    result: Dict[str, Any] = {}
    for pos, group in agg.groupby("Position"):
        entries: List[Dict[str, Any]] = []
        for _, row in group.iterrows():
            entries.append({
                "gene": row["Gene"],
                "transcripts": ','.join(row["Transcript"]),
                "protein_changes": ','.join(row["Protein_change"]),
                "ref_allele": ','.join(row["Ref_allele"])
            })
        result[pos] = entries[0] if len(entries) == 1 else entries

    return result


def get_read_bases_at_position(bam_path, chrom, pos) -> List[tuple]:
    """
    Return a list of tuples for each read covering a genomic position.

    Parameters
    ----------
    bam_path : str
        Path to an indexed BAM file.
    chrom : str
        Reference contig name (e.g. `chr13` or `13`).
    pos : int
        1-based genomic coordinate.

    Returns
    -------
    List[Tuple[str, str, int]]
        A list of `(read_id, base, q_score)` for each read that covers the
        specified position and has a base (i.e. not a deletion or ref-skip).
        `q_score` may be `None` if base qualities are not available.

    Notes
    -----
    This function uses `pysam.AlignmentFile.pileup` which operates on 0-based
    coordinates internally; the provided 1-based `pos` is converted accordingly.
    """
    results = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for pileupcolumn in bam.pileup(chrom, pos-1, pos, truncate=True):
            if pileupcolumn.pos == pos-1:
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        read = pileupread.alignment
                        base = read.query_sequence[pileupread.query_position]
                        q_score = read.query_qualities[pileupread.query_position]
                        results.append((read.query_name, base, q_score))
    return results


def encode_ref_alt_bases(base: str, ref: str) -> float:
    """
    Encode an observed read base relative to the reference allele.

    Parameters
    ----------
    base : str
        Single-character nucleotide observed in the read (for example: `A`, `C`, `G`, `T`).
    ref : str
        Single-character reference nucleotide at the genomic position.

    Returns
    -------
    str
        - `ref_allele` if `base == ref`.
        - `nonref_allele` if `base != ref`.
    """
    if base == ref:
        encoded = f"ref_allele"
    else:
        encoded = f"nonref_allele"
    return encoded


def write_reads_covering_positions(input_bam, output_bam, chrom, pos) -> None:
    """
    Write reads from `input_bam` that cover the genomic coordinate `chrom:pos` (1-based)
    into `output_bam`.

    Parameters
    ----------
    input_bam : str
        Path to the input BAM file.
    output_bam : str
        Path to the output BAM file to create.
    chrom : str
        Chromosome name (e.g. `chr13`).
    pos : int
        1-based genomic position to filter reads by.

    Notes
    -----
    The function performs a pileup at the specified position to collect read names
    that cover the base, then writes the
    full reads with those names to `output_bam`. The output BAM will be created
    using the input BAM as a template.
    """
    with pysam.AlignmentFile(input_bam, "rb") as bam:
        reads_pos = set()
        for pileupcolumn in bam.pileup(chrom, pos-1, pos, truncate=True):
            if pileupcolumn.pos == pos-1:
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        reads_pos.add(pileupread.alignment.query_name)
        bam.reset()
        with pysam.AlignmentFile(output_bam, "wb", template=bam) as out_bam:
            for read in bam.fetch(chrom):
                if read.query_name in reads_pos:
                    out_bam.write(read)


def analyze_hotspot_noise_for_coordinate(
    coord,
    samples_out,
    ref_samples,
    df_hspt,
    bam_dir,
    output_dir,
    hspt_info,
    bam_patterns=None,
    sample_id_col=None,
    bam_id_cols=None
):
    """
    Extract reads and compute sequencing noise statistics for a single coordinate.
    
    Parameters
    ----------
    coord : str
        Genomic coordinate, e.g. 'chr13:28018505'.
    samples_out : list of str
        Samples to exclude.
    ref_samples : pd.DataFrame
        DataFrame of reference samples.
    df_hspt : pd.DataFrame
        DataFrame containing hotspot metadata.
    bam_dir : str
        Directory containing BAM files.
    output_dir : str
        Directory to write output files.
    hspt_info : dict
        Dictionary of hotspots metadata to update with noise statistics.
    bam_patterns : list of str or callable, optional
        Templates for BAM filenames or callable for matching.
    sample_id_col : str or int, optional
        Column name or index in ref_samples for sample ID check.
    bam_id_cols : list of str or int, optional
        Column names or indices to construct BAM ID.

    Returns
    -------
    dict
        Dictionary of bases summaries for each processed sample.
    """
    if bam_patterns is None:
        bam_patterns = [
            "{bam_id}_T_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam",
            "{bam_id}_T.ont_reads.filtered.aligned.sorted.soft-clipped.bam"
        ]
    elif isinstance(bam_patterns, str):
        bam_patterns = [bam_patterns]

    chrom, pos = coord.split(":")[0], int(coord.split(":")[1])
    ref_allele = df_hspt[df_hspt['Position'] == coord]["Ref_allele"].drop_duplicates().values[0]
    print(f"Processing hotspot {coord}: {samples_out} will be excluded from the analysis.\n")
    
    coord_bases_table = {}
    last_bases_summary = None
    
    for ref_sample in ref_samples.itertuples():
        # Get the BAM ID
        if bam_id_cols is None:
            bam_id = '_'.join(str(val) for val in ref_sample[1:])
        else:
            vals = []
            for col in bam_id_cols:
                if isinstance(col, int):
                    vals.append(str(ref_sample[col + 1]))
                else:
                    vals.append(str(getattr(ref_sample, col)))
            bam_id = '_'.join(vals)
            
        # Get the Sample ID
        if sample_id_col is not None:
            if isinstance(sample_id_col, int):
                sample_id = ref_sample[sample_id_col + 1]
            else:
                sample_id = getattr(ref_sample, sample_id_col)
        else:
            sample_id = ref_sample[2] if len(ref_sample) > 2 else ref_sample[1]
            
        for bam in os.listdir(bam_dir):
            matched = False
            if callable(bam_patterns):
                matched = bam_patterns(bam, bam_id)
            else:
                for pattern in bam_patterns:
                    if "{bam_id}" in pattern:
                        if bam == pattern.format(bam_id=bam_id):
                            matched = True
                            break
                    else:
                        if bam == pattern:
                            matched = True
                            break
                            
            if not matched:
                continue
                
            print(ref_sample[1], sample_id)
            if sample_id in samples_out:
                print(f"Skipping sample {bam_id} for hotspot {coord} since the sample has a hotspot mutation.\n")
                continue
                
            sample_bam = os.path.join(bam_dir, bam)
            try:
                pysam.index(sample_bam)
            except Exception as e:
                print(f"Could not index {sample_bam}: {e}")
                continue
                
            bam_pos = os.path.join(output_dir, bam.replace(".bam", f".hotspot_{coord}.bam"))
            if not os.path.exists(bam_pos):
                print(f"Extracting reads covering {coord} from {sample_bam}.\nWriting them to {bam_pos}.\n")
                write_reads_covering_positions(sample_bam, bam_pos, chrom, pos)
                pysam.index(bam_pos)
                
            print(f"Extracting bases in reads in position {coord} from {bam_pos}.\n")
            hotspot_bases = get_read_bases_at_position(bam_pos, chrom, pos)
            bases_in_reads = pd.DataFrame(
                hotspot_bases,
                columns=["read_id", "base", "qscore"]
            )
            bases_in_reads.to_csv(os.path.join(output_dir, f"{bam_id}_{coord}_bases.csv"), index=False)
            
            bases_in_reads["base"] = bases_in_reads["base"].apply(lambda base: encode_ref_alt_bases(base, ref_allele))
            bases_summary = (
                bases_in_reads.groupby('base')['qscore']
                .agg(['count', 'mean'])
                .rename(columns={'count': 'read_count', 'mean': 'mean_qscore'})
            )
            
            if bases_summary is not None and not bases_summary.empty:
                bases_summary['read_count_ref_allele'] = bases_summary['read_count'].loc['ref_allele'] if 'ref_allele' in bases_summary.index else 0
                bases_summary['mean_qscore_ref_allele'] = bases_summary['mean_qscore'].loc['ref_allele'] if 'ref_allele' in bases_summary.index else 0
                bases_summary['read_count_nonref_allele'] = bases_summary['read_count'].loc['nonref_allele'] if 'nonref_allele' in bases_summary.index else 0
                bases_summary['mean_qscore_nonref_allele'] = bases_summary['mean_qscore'].loc['nonref_allele'] if 'nonref_allele' in bases_summary.index else 0
                
                # drop the row index 'ref_allele' if present
                if 'ref_allele' in bases_summary.index:
                    bases_summary.drop(index=['ref_allele'], inplace=True)
                bases_summary.drop(columns=['read_count', 'mean_qscore'], inplace=True)
                bases_summary.reset_index(inplace=True, drop=True)
                bases_summary['depth'] = bases_summary['read_count_ref_allele'] + bases_summary['read_count_nonref_allele']
                bases_summary['vaf'] = bases_summary['read_count_nonref_allele'] / bases_summary['depth']
                
                coord_bases_table[bam_id] = bases_summary
                last_bases_summary = bases_summary
                print(f"Summary stats of bases at {coord} for sample {bam_id}:\n{bases_summary}\r\n")
                
    print(f"Finished processing hotspot {coord}.\n")
    
    if coord_bases_table and last_bases_summary is not None:
        print([row.values.flatten() for row in coord_bases_table.values()])
        bases_table_df_coord = pd.DataFrame(
            index=coord_bases_table.keys(),
            columns=last_bases_summary.columns,
            data=[row.values.flatten() for row in coord_bases_table.values()]
        )
        bases_table_df_coord.to_csv(os.path.join(output_dir, f"hotspots_bases_summary_{coord}.csv"), index=True)
        print(f"Summary of bases at hotspot {coord} for all samples:\n{bases_table_df_coord}\r\n")
        
        # Compute and update noise statistics
        hspt_info[coord].update({
            "median_noise": bases_table_df_coord['vaf'].median(),
            "sd_noise": bases_table_df_coord['vaf'].std(),
            "threshold_5sd": bases_table_df_coord['vaf'].median() + 5 * bases_table_df_coord['vaf'].std(),
            "nb_samples": bases_table_df_coord.dropna(axis=0, how='any').shape[0]
        })
        print(f"Updated hotspot info for {coord} with noise statistics:\n{json.dumps(hspt_info[coord], indent=4)}\n")
        
        df_hspt_stats = pd.DataFrame.from_dict(hspt_info, orient="index")
        df_hspt_stats.index.rename("genomic_position_hg38", inplace=True)
        df_hspt_stats.to_csv(os.path.join(output_dir, f"hotspots_info_with_noise_stats.csv"), index=True)
        print(f"Finished processing hotspot {coord}.\n")
        
    return coord_bases_table


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract reads and compute sequencing noise statistics for hotspots.",
        fromfile_prefix_chars='@'
    )
    parser.add_argument("--bam_dir", help="Directory containing BAM files.")
    parser.add_argument("--hspt_csv", help="Path to the Hotspots CSV file.")
    parser.add_argument("--ref_panel", help="Path to the reference panel TSV file.")
    parser.add_argument("--bam_patterns", nargs="+", help="Templates for BAM filenames.")
    parser.add_argument("--sample_id_col", help="Column name or index in ref_samples for sample ID check.")
    parser.add_argument("--bam_id_cols", nargs="+", help="Column names or indices to construct BAM ID.")
    
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    bam_dir = args.bam_dir
    hspt_csv = args.hspt_csv
    ref_panel = args.ref_panel

    if not bam_dir or not hspt_csv or not ref_panel:
        raise ValueError(
            "bam_dir, hspt_csv, and ref_panel must be provided via command-line arguments or an @file."
        )

    ref_samples = pd.read_csv(ref_panel, sep="\t").drop(columns=["type", "caller"]).drop_duplicates()
    print(f"Reference samples used: {ref_samples}\n")

    df_hspt = pd.read_csv(hspt_csv, sep=",")
    exclude_samples = df_hspt.groupby("Position")["Sample"].apply(list).to_dict()
    print(f"Samples to exclude for each hotspot based on the presence of the hotspot mutation:\n{exclude_samples}\n")

    hspt_info = load_hotspots_metadata(hspt_csv)
    print(f"Information about each hotspot:\n{json.dumps(hspt_info, indent=4)}\n")

    bases_table = {}

    for coord, samples_out in exclude_samples.items():
        bases_table[coord] = analyze_hotspot_noise_for_coordinate(
            coord=coord,
            samples_out=samples_out,
            ref_samples=ref_samples,
            df_hspt=df_hspt,
            bam_dir=bam_dir,
            output_dir=os.path.dirname(ref_panel),
            hspt_info=hspt_info,
            bam_patterns=args.bam_patterns,
            sample_id_col=args.sample_id_col,
            bam_id_cols=args.bam_id_cols
        )
