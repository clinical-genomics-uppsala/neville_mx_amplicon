#import snakemake
from futils import *

df = FastxDataFrame(snakemake.input.fastq)
df.set_read_length
dshort = df.filter_too_short_reads()
write_frame_to_fastq(dshort, snakemake.output.fshort)
dlong = df.filter_too_long_reads()
write_frame_to_fastq(dlong, snakemake.output.flong)