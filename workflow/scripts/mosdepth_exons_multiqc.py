import gzip

with gzip.open(snakemake.input.bed, "rt") as fin, open(snakemake.output.tsv, "w") as fout:

    fout.write("# id: exon_coverage\n")
    fout.write("# section_name: Coverage per exon\n")
    fout.write("# description: Mean exon coverage from mosdepth\n")
    fout.write("# plot_type: table\n")
    fout.write("# pconfig:\n")
    fout.write("#   id: exon_coverage_table\n")
    fout.write("#   title: Coverage per exon\n")
    fout.write("#   namespace: Mosdepth\n")
    fout.write("Exon\tChromosome\tStart\tEnd\tMean coverage\n")

    for line in fin:
        chrom, start, end, exon, coverage = line.rstrip().split("\t")[:5]
        fout.write(
            f"{exon}\t{chrom}\t{start}\t{end}\t{coverage}\n"
        )
