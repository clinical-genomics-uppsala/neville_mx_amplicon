__author__ = "Camille Clouard"
__copyright__ = "Copyright 2024, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"

import os


rule varcall_clairs_to:
    input:
        bam="alignment/dorado_align/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam",
        bai="alignment/dorado_align/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam.bai",
        ref=config.get("ref_data"),
        bed=os.path.join(config.get("bed_files"), "amplicons.bed"),
    output:
        snv=temp("snv_indels/clairs_to/{sample}_{type}_snv.vcf.gz"),
        indel=temp("snv_indels/clairs_to/{sample}_{type}_indel.vcf.gz"),
    params:
        platform=config.get("clairs_to", {}).get("platform", ""),
        snv_min_af=config.get("clairs_to", {}).get("snv_min_af", 0.05),
        indel_min_af=config.get("clairs_to", {}).get("indel_min_af", 0.1),
        outdir=directory(lambda wildcards, output: os.path.dirname(output.snv)),
    resources:
        partition=config.get("clairs_to", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("clairs_to", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("clairs_to", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("clairs_to", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("clairs_to", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("clairs_to", {}).get("threads", config["default_resources"]["threads"])
    log:
        "snv_indels/clairs_to/{sample}_{type}_clairs_to.log",
    benchmark:
        repeat(
            "snv_indels/clairs_to/{sample}_{type}_clairs_to.benchmark.tsv",
            config.get("clairs_to", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("clairs_to", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Long-read somatic small variant calling in only tumor samples with ClairS-TO.
        """
    shell:
        """
        run_clairs_to --tumor_bam_fn {input.bam} --ref_fn {input.ref} --threads {resources.threads} --platform {params.platform} --output_dir {params.outdir} -s {wildcards.sample} --bed_fn {input.bed} --snv_min_af {params.snv_min_af} --indel_min_af {params.indel_min_af} --disable_verdict --snv_output_prefix {wildcards.sample}_{wildcards.type}_snv --indel_output_prefix {wildcards.sample}_{wildcards.type}_indel 2> {log}
        """


rule varcall_clairs_to_concat:
    input:
        snv="snv_indels/clairs_to/{sample}_{type}_snv.vcf.gz",
        indel="snv_indels/clairs_to/{sample}_{type}_indel.vcf.gz",
    output:
        all=temp(
            "snv_indels/clairs_to/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.clairs_to.vcf.gz"
        ),
    resources:
        partition=config.get("varcall_clairs_to_concat", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("varcall_clairs_to_concat", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("varcall_clairs_to_concat", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("varcall_clairs_to_concat", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("varcall_clairs_to_concat", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("varcall_clairs_to_concat", {}).get("threads", config["default_resources"]["threads"])
    container:
        config.get("varcall_clairs_to_concat", {}).get("container", config["default_container"])
    log:
        "snv_indels/clairs_to/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.clairs_to.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/clairs_to/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.clairs_to.benchmark.tsv",
            config.get("varcall_clairs_to_concat", {}).get("benchmark_repeats", 1),
        )
    message:
        """
        {rule}: Concatenate the output of ClairS-TO into a single VCF.
        """
    shell:
        """
        mkdir -p tmp/snv_indels/clairs_to
        bcftools concat -a -Oz -o tmp/{output.all} {input.snv} {input.indel} 2> {log}
        bcftools sort -Oz -o {output.all} tmp/{output.all}
        """


rule varcall_deepsomatic:
    input:
        bam="alignment/dorado_align/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam",
        bai= "alignment/dorado_align/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam.bai",
        ref=config.get("ref_data"),
        bed=config.get("deepsomatic", {}).get("bed_file", os.path.join(config.get("bed_files"),"amplicons.bed")),
    output:
        tmpdir=directory("snv_indels/deepsomatic/{sample}_{type}_tmp"),
        vcf="snv_indels/deepsomatic/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.deepsomatic.vcf.gz",
    params:
        sample=lambda wildcards: f"{wildcards.sample}",
        model="ONT_TUMOR_ONLY",
        extra=config.get("deepsomatic", {}).get("extra", ""),
        filter=config.get("deepsomatic", {}).get("filter", ""),
    log:
        "snv_indels/deepsomatic/{sample}_{type}_deepsomatic.log",
    benchmark:
        repeat(
            "snv_indels/deepsomatic/{sample}_{type}_deepsomatic.benchmark.tsv",
            config.get("deepsomatic", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("deepsomatic", {}).get("threads", config["default_resources"]["threads"])
    resources:
        partition=config.get("deepsomatic", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("deepsomatic", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("deepsomatic", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("deepsomatic", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepsomatic", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    container:
        config.get("deepsomatic", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Long-read somatic small variant calling in tumor samples only with DeepSomatic.
        """
    shell:
        """
        run_deepsomatic --model_type={params.model} --ref={input.ref} --reads_tumor={input.bam} --output_vcf={output.vcf} --sample_name_tumor={params.sample} --num_shards={resources.threads} --logging_dir={log} --intermediate_results_dir {output.tmpdir} --regions={input.bed} {params.extra} {params.filter}
        """


rule varcall_savana:
    input:
        bam="alignment/dorado_align/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam",
        bai="alignment/dorado_align/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam.bai",
        ref=config.get("ref_data"),
        bed=os.path.join(config.get("bed_files"), "amplicons.bed"),
    output:
        dummy="cnv_sv/savana/{sample}_{type}_savana.done",
        outdir=temp(directory("cnv_sv/savana/{sample}_{type}_savana_output")),
    params:
        sample=config.get("sample_id", "sample_T"),
        g1000_vcf="1000g_hg38",
        extra="",
        prefix=lambda w: f"{w.sample}_{w.type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.savana",
    log:
        "cnv_sv/savana/{sample}_{type}_savana.log",
    benchmark:
        repeat(
            "cnv_sv/savana/{sample}_{type}_savana.benchmark.tsv",
            config.get("savana", {}).get("benchmark_repeats", 1),
        )
    resources:
        partition=config.get("savana", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("savana", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("savana", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("savana", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("savana", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("savana", {}).get("threads", config["default_resources"]["threads"])
    container:
        config.get("savana", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Long-read somatic small variant and SV calling in tumor samples only with SAVANA.
        """
    shell:
        """
        savana --help
        savana to --tumour {input.bam} --outdir {output.outdir} --ref {input.ref} --g1000_vcf {params.g1000_vcf} --chromosomes 2 5 13 15 17 &> {log}
        touch {output.dummy}
        """