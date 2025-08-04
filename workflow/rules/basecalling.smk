__author__ = "Camille Clouard"
__copyright__ = "Copyright 2024, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"

import os

rule basecalling_dorado:
    input:
        pod5 = os.path.join(config.get("runfolder"),
            config.get("raw_data")
        ),
        ref_data = config.get("ref_data")
    output:
        bam = temp("basecalling/dorado/{sample}_{type}_reads.ont_adapt_trim.bam")
    params:
        dir_models = config.get("dir_models"),  # os.path.abspath(os.pardir),
        dorado_model = config.get("dorado_model"),
        dorado_options = "--device cuda:all --trim adapters"  # --trim all
    resources:
        partition=config.get("dorado_basecalling",{}).get("partition",config["default_resources"]["partition"]),
        time=config.get("dorado_basecalling",{}).get("time",config["default_resources"]["time"]),
        gres=config.get("dorado_basecalling",{}).get("gres"),
        threads=config.get("dorado_basecalling",{}).get("threads",config["default_resources"]["threads"]),
        mem_mb=config.get("dorado_basecalling",{}).get("mem_mb",config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("dorado_basecalling",{}).get("mem_per_cpu",config["default_resources"]["mem_per_cpu"]),
        slurm_extra=config.get("dorado_basecalling",{}).get("slurm_extra"),
    threads: config.get("dorado_basecalling", {}).get("threads", config["default_resources"]["threads"]),
    container:
        config.get("dorado", {}).get("container", config["default_container"])
    log:
        "basecalling/dorado/{sample}_{type}_reads.ont_adapt_trim.bam.log"
    benchmark:
        repeat(
            "basecalling/dorado/{sample}_{type}_reads.ont_adapt_trim.bam.benchmark.tsv",
            config.get("dorado_basecalling", {}).get("benchmark_repeats", 1)
        )
    message:
        "{rule}: Basecalling with dorado from POD5 files. ONT adapters will be trimmed."
    shell:
        """
        echo "Dorado executed from $( which dorado )"
        
        echo "Executing dorado basecalling in {input.pod5} with options '{params.dorado_options}'"
        echo "and model {params.dorado_model}"
        
        dorado basecaller {params.dir_models}/{params.dorado_model} {params.dorado_options} {input.pod5}/ > {output.bam} 2>> {log}
        """

rule duplex_basecalling_dorado:
    input:
        pod5 = os.path.join(config.get("runfolder"),
            config.get("raw_data")
        ),
        ref_data = config.get("ref_data")
    output:
        bam = temp("basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.bam")
    params:
        dir_models=config.get("dir_models"),
        dorado_model=config.get("dorado_model"),
        dorado_options="--device cuda:all",
        trim_options="--no-trim-primers --sequencing-kit SQK-LSK114"
    resources:
        partition=config.get("dorado_basecalling",{}).get("partition",config["default_resources"]["partition"]),
        time=config.get("dorado_basecalling",{}).get("time",config["default_resources"]["time"]),
        gres=config.get("dorado_basecalling",{}).get("gres"),
        threads=config.get("dorado_basecalling",{}).get("threads",config["default_resources"]["threads"]),
        mem_mb=config.get("dorado_basecalling",{}).get("mem_mb",config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("dorado_basecalling",{}).get("mem_per_cpu",config["default_resources"]["mem_per_cpu"]),
        slurm_extra=config.get("dorado_basecalling",{}).get("slurm_extra"),
    threads: config.get("dorado_basecalling", {}).get("threads", config["default_resources"]["threads"]),
    benchmark:
        repeat(
            "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.bam.benchmark.tsv",
            config.get("dorado_basecalling", {}).get("benchmark_repeats", 1)
        )
    container:
        config.get("dorado",{}).get("container",config["default_container"])
    log:
        "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.bam.log"
    message:
        "{rule}: Duplex basecalling with dorado from POD5 files. ONT adapters will be trimmed."
    shell:
        """
        echo "Dorado executed from $( which dorado )"

        echo "Executing dorado duplex basecalling in {input.pod5} with options '{params.dorado_options}'"
        echo "and model {params.dorado_model}"

        dorado duplex {params.dir_models}/{params.dorado_model} {params.dorado_options} {input.pod5}/ | dorado trim {params.trim_options} > {output.bam} 2> {log}
        """

rule basecalling_bam2fastq:
    input:
        bam = "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.bam",
    output:
        fastq = temp("basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq")
    resources:
        partition=config.get("samtools",{}).get("partition",config["default_resources"]["partition"]),
        time=config.get("samtools",{}).get("time",config["default_resources"]["time"]),
        threads=config.get("samtools",{}).get("threads",config["default_resources"]["threads"]),
        mem_mb=config.get("samtools",{}).get("mem_mb",config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("samtools",{}).get("mem_per_cpu",config["default_resources"]["mem_per_cpu"]),
    threads: config.get("samtools", {}).get("threads", config["default_resources"]["threads"]),
    benchmark:
        repeat(
            "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq.benchmark.tsv",
            config.get("samtools", {}).get("benchmark_repeats", 1)
        )
    log:
        "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq.log"
    container:
        config.get("samtools", {}).get("container", config["default_container"])
    message:
        "{rule}: Convert unaligned BAM file to FASTQ with samtools"
    shell:
        """
        samtools fastq {input.bam} > {output.fastq} 2> {log}
        """

rule basecalling_compress_fastq:
    input:
        fastq = "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq"
    output:
        fastqgz = temp("basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq.gz")
    resources:
        partition=config.get("samtools",{}).get("partition",config["default_resources"]["partition"]),
        time=config.get("samtools",{}).get("time",config["default_resources"]["time"]),
        threads=config.get("samtools",{}).get("threads",config["default_resources"]["threads"]),
        mem_mb=config.get("samtools",{}).get("mem_mb",config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("samtools",{}).get("mem_per_cpu",config["default_resources"]["mem_per_cpu"]),
    threads: config.get("samtools", {}).get("threads", config["default_resources"]["threads"]),
    benchmark:
        repeat(
            "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq.gz.benchmark.tsv",
            config.get("samtools", {}).get("benchmark_repeats", 1)
        )
    log:
        "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq.gz.log"
    container:
        config.get("basecalling_compress_fastq", {}).get("container", config["default_container"])
    message:
        "{rule}: Compress FASTQ file with basecalled data."
    shell:
        """
        rm -f {output.fastqgz}
        gzip -f {input.fastq} > {output.fastqgz} 2> {log}
        """