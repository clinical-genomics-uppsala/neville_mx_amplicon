from xopen import xopen
import pysam
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


if __name__ == "__main__":
    fastx_path = "/pipeline_pool_amplicon/run/"\
                 "20240923_1256_MN45214_ASA641_69751c45/data/reads.ont_adapt_trim.filtered.out.fastq.gz"
    fastqshort = "/home/camille/ampliconthemato/pipeline_pool_amplicon/run/"\
                 "20240923_1256_MN45214_ASA641_69751c45/data/reads.ont_adapt_trim.filtered.out.short.fastq.gz"
    df = FastxDataFrame(fastx_path)
    df.set_read_length
    df.set_read_phred_score
    dreads = df.get_frame()
    print(dreads["phred_scores"])
    print(dreads["length"])
    print(dreads["q_score"])
    dshort = df.filter_too_short_reads()
    dlong = df.filter_too_long_reads()
    print(dshort.head()["length"])
    print(dlong.head()["length"])
    write_frame_to_fastq(dshort, fastqshort)
