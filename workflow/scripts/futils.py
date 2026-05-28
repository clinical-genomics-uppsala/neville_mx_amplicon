from xopen import xopen
import pysam

class FastxDataFrame(object):
    """
    A lightweight, streaming FASTQ parser and splitter.
    Replaces the memory-heavy Pandas DataFrame implementation to prevent OOM errors
    and process large FASTQ files efficiently in a single pass.
    """
    def __init__(self, fastx_path):
        self.fastx_path = fastx_path

    def split_by_length(self, short_path, long_path, min_length=2000, max_length=4000):
        """
        Streams the input FASTQ file and writes reads directly to short/long output files
        based on length thresholds in a single efficient pass.
        """
        with pysam.FastxFile(self.fastx_path) as infile, \
             xopen(short_path, "w") as out_short, \
             xopen(long_path, "w") as out_long:
            for read in infile:
                rlen = len(read.sequence)
                # Format quality and comments if present
                comment_part = f" {read.comment}" if read.comment else ""
                fastq_record = f"@{read.name}{comment_part}\n{read.sequence}\n+\n{read.quality}\n"
                
                if rlen < min_length:
                    out_short.write(fastq_record)
                elif rlen > max_length:
                    out_long.write(fastq_record)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Split FASTQ file into too-short and too-long reads.")
    parser.add_argument("input_fastq", help="Path to input FASTQ file")
    parser.add_argument("output_short", help="Path to output too-short FASTQ file")
    parser.add_argument("output_long", help="Path to output too-long FASTQ file")
    parser.add_argument("--min-len", type=int, default=2000, help="Min length threshold")
    parser.add_argument("--max-len", type=int, default=4000, help="Max length threshold")
    args = parser.parse_args()

    df = FastxDataFrame(args.input_fastq)
    df.split_by_length(args.output_short, args.output_long, args.min_len, args.max_len)
