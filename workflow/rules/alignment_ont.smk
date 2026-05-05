__author__ = "Camille Clouard"
__copyright__ = "Copyright 2026, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"

import os
from snakemake.logging import logger

logger.info(f"\n{workflow.snakefile} is being parsed")


rule alignment_ont_dorado_align:
    input:
        fastqgz="prealignment/filtlong/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.fastq.gz",
        ref_data=config.get("ref_data"),
    output:
        bam=temp("alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.bam"),
    params:
        runid=config.get("runid", ""),
        extra=config.get("alignment_ont_dorado_align", {}).get("extra", ""),
    resources:
        partition=config.get("alignment_ont_dorado_align", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("alignment_ont_dorado_align", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("alignment_ont_dorado_align", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("alignment_ont_dorado_align", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("alignment_ont_dorado_align", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
    threads: config.get("alignment_ont_dorado_align", {}).get("threads", config["default_resources"]["threads"])
    container:
        config.get("alignment_ont_dorado_align", {}).get("container", config["default_container"])
    log:
        "alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.log",
    benchmark:
        repeat(
            "alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.benchmark.tsv",
            config.get("alignment_ont_dorado_align", {}).get("benchmark_repeats", 1),
        )
    message:
        "{rule}: Align reads with dorado and minimap2"
    shell:
        """
        echo "Dorado executed from $( which dorado )"

        echo "Executing dorado aligning of {input.fastqgz} with reference genome '{input.ref_data}'"

        dorado aligner {params.extra} {input.ref_data} {input.fastqgz} > {output.bam} 2> {log}

        if [ {params.runid} == "ABC123" ]; then
            echo "\nIntegration test detected inside dorado_align rule, overwriting the output BAM file with a properly aligned BAM file" >> {log}
            samtools view -b test_data/preprocessed/Mtest_D25-test007_T_reads.ont_adapt_trim.filtered.aligned.bam > {output.bam} 2>> {log}
        fi
        """


rule alignment_ont_bam_sort:
    input:
        "alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.bam",
    output:
        temp("alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.bam"),
    resources:
        partition=config.get("alignment_ont_bam_sort", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("alignment_ont_bam_sort", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("alignment_ont_bam_sort", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("alignment_ont_bam_sort", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("alignment_ont_bam_sort", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("alignment_ont_bam_sort", {}).get("threads", config["default_resources"]["threads"])
    container:
        config.get("alignment_ont_bam_sort", {}).get("container", config["default_container"])
    log:
        "alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.log",
    benchmark:
        repeat(
            "alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.benchmark.tsv",
            config.get("alignment_ont_bam_sort", {}).get("benchmark_repeats", 1),
        )
    message:
        "{rule}: Sort aligned reads with samtools"
    wrapper:
        "0.2.0/bio/samtools/sort"


rule alignment_ont_bam_index:
    input:
        "alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.bam",
    output:
        temp("alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.bam.bai"),
    resources:
        partition=config.get("alignment_ont_bam_index", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("alignment_ont_bam_index", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("alignment_ont_bam_index", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("alignment_ont_bam_index", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("alignment_ont_bam_index", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("alignment_ont_bam_index", {}).get("threads", config["default_resources"]["threads"])
    container:
        config.get("alignment_ont_bam_index", {}).get("container", config["default_container"])
    log:
        "alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.bam.bai.log",
    benchmark:
        repeat(
            "alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.bam.bai.benchmark.tsv",
            config.get("alignment_ont_bam_index", {}).get("benchmark_repeats", 1),
        )
    message:
        "{rule}: Index the aligned and sorted reads with samtools"
    wrapper:
        "0.2.0/bio/samtools/index"


rule alignment_ont_bam_softclip:
    input:
        bam="alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.bam",
        amplibed=os.path.join(config.get("bed_files"), "amplicons-primers.bed"),
    output:
        bamclip=temp(
            "alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam"
        ),
        bamclipidx=temp(
            "alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam.bai"
        ),
        summary=temp("results/qc/samtools_ampliconclip/{experiment}_{sample}_{type}_ampliconclip.txt"),
    resources:
        partition=config.get("alignment_ont_bam_softclip", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("alignment_ont_bam_softclip", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("alignment_ont_bam_softclip", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("alignment_ont_bam_softclip", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("alignment_ont_bam_softclip", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("alignment_ont_bam_softclip", {}).get("threads", config["default_resources"]["threads"])
    log:
        "alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.log",
    benchmark:
        repeat(
            "alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.benchmark.tsv",
            config.get("alignment_ont_bam_softclip", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("alignment_ont_bam_softclip", {}).get("container", config["default_container"])
    message:
        "{rule}: Mark the primers as soft-clipped bases with samtools"
    shell:
        """
        samtools ampliconclip -o {output.bamclip} -b {input.amplibed} -f {output.summary} --soft-clip --both-ends --clipped {input.bam} 2> {log}
        samtools index {output.bamclip} 2>> {log}
        """


rule alignment_ont_split_bam_by_target:
    input:
        bam="alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam",
        bai="alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam.bai",
        amplibed=os.path.join(config.get("bed_files"), "{target}.bed"),
    output:
        bam=temp(
            "alignment/dorado_align/{experiment}_{sample}_{type}_{target}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam"
        ),
        bai=temp(
            "alignment/dorado_align/{experiment}_{sample}_{type}_{target}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam.bai"
        ),
        txt="alignment/dorado_align/{experiment}_{sample}_{type}_{target}_reads.txt",
    resources:
        partition=config.get("samtools", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("samtools", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("samtools", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("samtools", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("samtools", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("samtools", {}).get("threads", config["default_resources"]["threads"])
    benchmark:
        repeat(
            "alignment/dorado_align/{experiment}_{sample}_{type}_{target}_split.bam.benchmark.tsv",
            config.get("samtools", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("samtools", {}).get("container", config["default_container"])
    log:
        "alignment/dorado_align/{experiment}_{sample}_{type}_{target}_split.bam.log",
    message:
        "{rule}: Write one BAM file per target"
    shell:
        """
        samtools view -L {input.amplibed} -bo {output.bam} {input.bam} 2> {log}
        samtools index {output.bam} 2>> {log}
        samtools view {output.bam} | cut -d$'\t' -f1 > {output.txt} 2>> {log}
        """


rule alignment_ont_split_bam_target_j3:
    input:
        txt=expand(
            "alignment/dorado_align/{{experiment}}_{{sample}}_{{type}}_{target}_reads.txt",
            target=config.get("amplicons") + config.get("extra_regions"),
        ),
    output:
        bam=temp(
            "alignment/dorado_align/{experiment}_{sample}_{type}_TP53_J3_only_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam"
        ),
        bai=temp(
            "alignment/dorado_align/{experiment}_{sample}_{type}_TP53_J3_only_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam.bai"
        ),
        txt="alignment/dorado_align/{experiment}_{sample}_{type}_TP53_J3_only_reads.txt",
    params:
        inbam=lambda wildcards, input: os.path.join(
            os.path.dirname(input[0]),
            f"{wildcards.experiment}_{wildcards.sample}_{wildcards.type}_TP53_D2+J3_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam",
        ),
        intxt=lambda wildcards, input: os.path.join(
            os.path.dirname(input[0]), f"{wildcards.experiment}_{wildcards.sample}_{wildcards.type}_TP53_D2_only_reads.txt"
        ),
    resources:
        partition=config.get("alignment_ont_split_bam_target_j3", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("alignment_ont_split_bam_target_j3", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("alignment_ont_split_bam_target_j3", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("alignment_ont_split_bam_target_j3", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("alignment_ont_split_bam_target_j3", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("alignment_ont_split_bam_target_j3", {}).get("threads", config["default_resources"]["threads"])
    benchmark:
        repeat(
            "alignment/dorado_align/{experiment}_{sample}_{type}_TP53_J3_only_split.bam.benchmark.tsv",
            config.get("alignment_ont_split_bam_target_j3", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("alignment_ont_split_bam_target_j3", {}).get("container", config["default_container"])
    log:
        "alignment/dorado_align/{experiment}_{sample}_{type}_TP53_J3_only_split.bam.log",
    message:
        "{rule}: Handle the specific case of totally overlapping amplicons D2+J3 and write one BAM file for J3"
    shell:
        """
        samtools view -bo {output.bam} -N ^{params.intxt} {params.inbam} 2> {log}
        samtools index {output.bam} 2>> {log}
        samtools view {output.bam} | cut -d$'\t' -f1 > {output.txt} 2>> {log}
        """
