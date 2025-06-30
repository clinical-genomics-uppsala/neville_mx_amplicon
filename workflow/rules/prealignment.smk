__author__ = "Camille Clouard"
__copyright__ = "Copyright 2024, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"

rule filtlong:
    input:
        "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq.gz"
    output:
        fastq=temp("prealignment/filtlong/{sample}_{type}_reads.ont_adapt_trim.filtered.fastq.gz"),
        # info="{experiment_name}/{sample_id}/{run_date}_{run_id}/results/filtlong.txt"
    params:
        min_length = config.get("filtlong", {}).get("min_length", 150),
        max_length = config.get("filtlong", {}).get("max_length", 8000)
    log:
        "prealignment/filtlong/{sample}_{type}.filtlong.log"
    resources:
        partition=config.get("filtlong",{}).get("partition",config["default_resources"]["partition"]),
        time=config.get("filtlong",{}).get("time",config["default_resources"]["time"]),
        threads=config.get("filtlong",{}).get("threads",config["default_resources"]["threads"]),
        mem_mb=config.get("filtlong",{}).get("mem_mb",config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("filtlong",{}).get("mem_per_cpu",config["default_resources"]["mem_per_cpu"]),
    container:
        config.get("filtlong", {}).get("container", config["default_container"])
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
        fastq1="prealignment/filtlong/{sample}_{type}_reads.ont_adapt_trim.filtered.fastq.gz",
        fastq2="basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq.gz"
    output:
        txt=temp("prealignment/filtlong/{sample}_{type}_filtered_read_id.txt"),
        bam=temp("prealignment/filtlong/{sample}_{type}_reads.ont_adapt_trim.filtered.out.bam"),
        fastq=temp("prealignment/filtlong/{sample}_{type}_reads.ont_adapt_trim.filtered.out.fastq.gz")
    resources:
        partition = config.get("filtlong",{}).get("partition",config["default_resources"]["partition"]),
        time = config.get("filtlong",{}).get("time",config["default_resources"]["time"]),
        threads = config.get("filtlong",{}).get("threads",config["default_resources"]["threads"]),
        mem_mb = config.get("filtlong",{}).get("mem_mb",config["default_resources"]["mem_mb"]),
        mem_per_cpu = config.get("filtlong",{}).get("mem_per_cpu",config["default_resources"]["mem_per_cpu"])
    container:
        config.get("samtools", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Fetch the reads that are filtered out with Filtlong.
        """
    shell:
        """
        samtools view {input.fastq1} | cut -d$'\t' -f1 > {output.txt}
        samtools view -bo {output.bam} -N ^{output.txt} {input.fastq2}
        samtools fastq {output.bam} | gzip > {output.fastq}
        """

rule fetch_too_short_too_long_reads:
    input:
        fastq = "prealignment/filtlong/{sample}_{type}_reads.ont_adapt_trim.filtered.out.fastq.gz"
    output:
        fshort = temp("prealignment/filtlong/{sample}_{type}_reads.ont_adapt_trim.filtered.out.short.fastq.gz"),
        flong = temp("prealignment/filtlong/{sample}_{type}_reads.ont_adapt_trim.filtered.out.long.fastq.gz"),

    resources:
        partition=config["default_resources"]["partition"],
        time=config["default_resources"]["time"],
        threads=config["default_resources"]["threads"],
        mem_mb=config["default_resources"]["mem_mb"],
        mem_per_cpu=config["default_resources"]["mem_per_cpu"]
    container:
        config["default_container"]
    message:
        """
        {rule}: Fetch the reads that are too short or too long and write them to compressed FASTQ format.
        """
    script:
        "../scripts/fetch_too_short_too_long_reads.py"
