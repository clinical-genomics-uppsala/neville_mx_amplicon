from futils import FastxDataFrame

df = FastxDataFrame(snakemake.input.fastq)
df.split_by_length(snakemake.output.fshort, snakemake.output.flong)
