__author__ = "Camille Clouard"
__copyright__ = "Copyright 2024, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"

from snakemake.logging import logger

logger.info(f"\n{workflow.snakefile} is being parsed")


def get_fastqfilter_input(wildcards):
    if config.get("data_type", "pod5") == "fastq":
        return config.get("input_fastq").format(
            experiment=wildcards.experiment,
            sample=wildcards.sample,
            type=wildcards.type
        )
    return "basecalling/dorado_duplex/{experiment}_{sample}_{type}_reads.ont_adapt_trim.fastq.gz".format(
        experiment=wildcards.experiment,
        sample=wildcards.sample,
        type=wildcards.type
    )

def get_nanoparse_input(wildcards):
    if config.get("nanoMonitor", {}).get("use_fastqfilter", False):
        return "prealignment/nanomonitor/{experiment}_{sample}_{type}_reads.fastqfilter.fastq.gz".format(
            experiment=wildcards.experiment,
            sample=wildcards.sample,
            type=wildcards.type
        )
    return get_fastqfilter_input(wildcards)

def get_filtlong_input_chained(wildcards):
    if config.get("nanoMonitor", {}).get("use_nanoparse", False):
        return "prealignment/nanomonitor/{experiment}_{sample}_{type}_reads.nanoparse.fastq.gz".format(
            experiment=wildcards.experiment,
            sample=wildcards.sample,
            type=wildcards.type
        )
    return get_nanoparse_input(wildcards)


if config.get("nanoMonitor", {}).get("use_fastqfilter", False):
    rule nanomonitor_fastqfilter:
        input:
            fastq=get_fastqfilter_input,
        output:
            fastq=temp("prealignment/nanomonitor/{experiment}_{sample}_{type}_reads.fastqfilter.fastq.gz"),
        params:
            bin=config.get("nanoMonitor", {}).get("fastqfilter_bin", "bin/fastqfilter"),
            extra=config.get("nanoMonitor", {}).get("fastqfilter_extra", ""),
        shell:
            "{params.bin} -i {input.fastq} -o {output.fastq} {params.extra}"


if config.get("nanoMonitor", {}).get("use_nanoparse", False):
    rule nanomonitor_nanoparse:
        input:
            fastq=get_nanoparse_input,
        output:
            fastq=temp("prealignment/nanomonitor/{experiment}_{sample}_{type}_reads.nanoparse.fastq.gz"),
        params:
            bin=config.get("nanoMonitor", {}).get("nanoparse_bin", "bin/nanoparse"),
            extra=config.get("nanoMonitor", {}).get("nanoparse_extra", ""),
        shell:
            "{params.bin} -i {input.fastq} -o {output.fastq} {params.extra}"


rule filtlong:
    input:
        get_filtlong_input_chained,
    output:
        fastq=temp("prealignment/filtlong/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.fastq.gz"),
    params:
        min_length=config.get("filtlong", {}).get("min_length", 150),
        max_length=config.get("filtlong", {}).get("max_length", 8000),
    log:
        "prealignment/filtlong/{experiment}_{sample}_{type}.filtlong.log",
    resources:
        partition=config.get("filtlong", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("filtlong", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("filtlong", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("filtlong", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("filtlong", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("filtlong", {}).get("threads", config["default_resources"]["threads"])
    container:
        config.get("filtlong", {}).get("container", config["default_container"])
    benchmark:
        repeat(
            "prealignment/filtlong/{experiment}_{sample}_{type}.filtlong.benchmark.tsv",
            config.get("filtlong", {}).get("benchmark_repeats", 1),
        )
    message:
        """
        {rule}: Filter reads with Filtlong. Keep reads of good quality
        that have a length >= {params.min_length} and <= {params.max_length}.
        """
    shell:
        """
        filtlong --min_length {params.min_length} --max_length {params.max_length} --min_mean_q 30 --keep_percent 75 {input}  2> {log} | gzip > {output.fastq}
        """


rule fetch_filtered_reads:
    input:
        fastq1="prealignment/filtlong/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.fastq.gz",
        fastq2="basecalling/dorado_duplex/{experiment}_{sample}_{type}_reads.ont_adapt_trim.fastq.gz",
    output:
        txt=temp("prealignment/filtlong/{experiment}_{sample}_{type}_filtered_read_id.txt"),
        bam=temp("prealignment/filtlong/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.out.bam"),
        fastq=temp("prealignment/filtlong/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.out.fastq.gz"),
    resources:
        partition=config.get("filtlong", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("filtlong", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("filtlong", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("filtlong", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("filtlong", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("filtlong", {}).get("threads", config["default_resources"]["threads"])
    log:
        "prealignment/filtlong/{experiment}_{sample}_{type}.fetch_filtered_reads.log",
    benchmark:
        repeat(
            "prealignment/filtlong/{experiment}_{sample}_{type}.fetch_filtered_reads.benchmark.tsv",
            config.get("filtlong", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("samtools", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Fetch the reads that are filtered out with Filtlong.
        """
    shell:
        """
        samtools view {input.fastq1} | cut -d$'\t' -f1 > {output.txt}
        samtools view -bo {output.bam} -N ^{output.txt} {input.fastq2} 2> {log}
        samtools fastq {output.bam} | gzip > {output.fastq} 2>> {log}
        """


rule fetch_too_short_too_long_reads:
    input:
        fastq="prealignment/filtlong/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.out.fastq.gz",
    output:
        fshort=temp("prealignment/filtlong/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.out.short.fastq.gz"),
        flong=temp("prealignment/filtlong/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.out.long.fastq.gz"),
    resources:
        partition=config["default_resources"]["partition"],
        time=config["default_resources"]["time"],
        threads=config["default_resources"]["threads"],
        mem_mb=config["default_resources"]["mem_mb"],
        mem_per_cpu=config["default_resources"]["mem_per_cpu"],
    threads: config["default_resources"]["threads"]
    log:
        "prealignment/filtlong/{experiment}_{sample}_{type}.fetch_too_short_too_long_reads.log",
    benchmark:
        repeat(
            "prealignment/filtlong/{experiment}_{sample}_{type}.fetch_too_short_too_long_reads.benchmark.tsv",
            config.get("filtlong", {}).get("benchmark_repeats", 1),
        )
    container:
        config["default_container"]
    message:
        """
        {rule}: Fetch the reads that are too short or too long and write them to compressed FASTQ format.
        """
    script:
        "../scripts/fetch_too_short_too_long_reads.py"
