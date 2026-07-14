from xopen import xopen
import pysam
import numpy as np
import pandas as pd


class FastxDataFrame(object):
    """
    Class to read a FASTQ file and convert it into a pandas DataFrame.
    :param fastx_path: path to the FASTQ file
    :return: a pandas DataFrame containing information from the reads in the following columns:
    - sequence: the nucleotide sequence of the read
    - comment: the comment associated with the read (if any)
    - quality: the quality string of the read
    - phred_scores: the quality scores of the read as a list of integers
    - length: the length of the read (calculated property)
    - q_score: the average quality score of the read (calculated property)
    The class provides methods to filter reads based on their length and to write the DataFrame back to a FASTQ file.
    """
    def __init__(self, fastx_path):
        self.fastxf = pysam.FastxFile(fastx_path)
        fastxdic = dict()
        for read in self.fastxf:
            fastxdic[read.name] = [read.sequence,
                                   read.comment,
                                   read.quality,
                                   read.get_quality_array(offset=33)]
        self.fastxdf = pd.DataFrame.from_dict(fastxdic, orient="index",
                                              columns=["sequence",
                                                       "comment",
                                                       "quality",
                                                       "phred_scores"],
                                              dtype=object
                                              )

    def get_frame(self):
        return self.fastxdf

    @property
    def set_read_length(self):
        self.fastxdf = self.fastxdf.assign(length=lambda x: x.sequence.str.len())

    @property
    def set_read_phred_score(self):  # use pysam.FastqProxy instead?
        """
        https://nanoporetech.com/support/software/data-analysis/where-can-i-find-out-more-about-quality-scores:
        Per-base quality scores are stored together with the base sequence in FASTQ files output
        by the basecalling algorithms and are then encoded in the Sanger format using ASCII characters
        with values of 33 to 126 (up to 93 ASCII character values).
        """
        self.fastxdf["q_score"] = self.fastxdf["phred_scores"].map(lambda x: round(sum(x)/len(x)))

    def filter_too_short_reads(self, min_length=2000):
        if "length" not in self.fastxdf.columns:
            self.set_read_length
        return self.fastxdf[(self.fastxdf["length"] < min_length)]

    def filter_too_long_reads(self, max_length=4000):
        if "length" not in self.fastxdf.columns:
            self.set_read_length
        return self.fastxdf[(self.fastxdf["length"] > max_length)]


def write_frame_to_fastq(dframe, pathout):
    with xopen(pathout, "w") as fout:
        for row in dframe.itertuples():
            fout.write(f"@{row.Index}\n{row.sequence}\n+\n{row.quality}\n")


# Function than extracts the average basecalling quality per base in a limited region from a BAM file using pysam
# and returns a tsv file with chrom and pos columns, followed by a column with the average quality per base in the region
def extract_average_quality_from_bam(bam_path, bed_path, tsv_out_path):
    """
    Extracts the average basecalling quality per base in a limited region from a BAM file using pysam
    and returns a tsv file with the usual BED columns followed by a column with the average quality per base in the region.
    :param bam_path: path to the BAM file
    :param bed_path: path to the BED file containing the regions of interest
    :param tsv_out_path: path to the output TSV file
    """
    # Read the BED file into a pandas DataFrame
    df_bed = pd.read_csv(bed_path, sep='\t', header=None, names=["chrom", "start", "end", "name", "pool"])

    # Open the BAM file using pysam
    bamfile = pysam.AlignmentFile(bam_path, "rb")

    # Prepare a list to hold the results
    results = []

    # Iterate over each region in the BED file
    for target in df_bed.itertuples():
        target_name = target.name
        target_start = target.start
        target_end = target.end
        target_chrom = target.chrom

        dict_pos_qscore = dict([(int(pos), list()) for pos in np.arange(start=target_start, stop=target_end, step=1)])

        # Fetch reads overlapping the region of interest
        reads = bamfile.fetch(target_chrom, target_start, target_end)

        # Initialize a list to hold quality scores for each base in the region
        quality_scores = []

        # Iterate over each read and extract quality scores for bases within the region
        for read in reads:
            print(read.to_string())
            if not read.is_unmapped:
                print(read.query_length)
                print(read.get_reference_positions())
                print(len(read.get_reference_positions()))
                for i in range(read.query_length):
                    try:
                        ref_pos = read.get_reference_positions()[i]
                        if ref_pos in dict_pos_qscore.keys():
                            dict_pos_qscore[int(ref_pos)].append(read.query_qualities[i])
                        else:
                            dict_pos_qscore[int(ref_pos)] = [read.query_qualities[i]]
                    except IndexError:
                        continue
                print(dict_pos_qscore)
                break

        # Calculate average quality score for the region
        dict_pos_avg_quality = dict()
        for pos, qscores in dict_pos_qscore.items():
            dict_pos_avg_quality[pos] = round(sum(qscores) / len(qscores), 2) if qscores != [] else 0.0

    # Convert results to a DataFrame and save as TSV
    df_results = pd.DataFrame.from_dict(dict_pos_avg_quality, orient='index', columns=["avg_base_quality"])
    #df_results.index = [f"{target_chrom}:{pos}" for pos in df_results.index]
    df_results["Chr"] = target_chrom.strip("chr")
    df_results["Pos"] = df_results.index
    df_results.reset_index()
    print(df_results)
    df_results.to_csv(tsv_out_path, sep='\t', index=False)

    # Close the BAM file
    bamfile.close()

if __name__ == "__main__":
    pass
    bam_path = "/home/camille/ampliconthemato/neville_mx_amplicon/hotspot_noise_analysis/M11_D23-09080_T.ont_reads.filtered.aligned.sorted.soft-clipped.bam"
    bed_path = "/home/camille/ampliconthemato/amplic_ont_hemato/data/primer_data/FLT3_ITD_3kb_1.bed"
    tsv_out = "/home/camille/ampliconthemato/neville_mx_amplicon/hotspot_noise_analysis/M11_D23-09080_T.FLT3-ITD.qscores.tsv"

    extract_average_quality_from_bam(bam_path, bed_path, tsv_out)
