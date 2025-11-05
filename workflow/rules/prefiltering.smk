__author__ = "Camille Clouard"
__copyright__ = "Copyright 2024, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"

from snakemake.logging import logger

logger.info(f"\n{workflow.snakefile} is being parsed")


rule filtering_bcftools_view:
    input:
        vcf="snv_indels/{caller}/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.{caller}.vcf.gz",
        tbi="snv_indels/{caller}/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.{caller}.vcf.gz.tbi",
    output:
        vcf=temp("snv_indels/{caller}/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.{caller}.bcftools_view.vcf.gz"),
    params:
        extra=lambda wildcards: config.get("bcftools_view",{}) \
                                    .get("filter", {}) \
                                    .get(f"{wildcards.caller}",{}).get("extra", "")
    log:
        "snv_indels/{caller}/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.{caller}.bcftools_view.vcf.log",
    benchmark:
        repeat("snv_indels/{caller}/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.{caller}.bcftools_view.vcf.benchmark.tsv",
            config.get("bcftools_view",{}).get("benchmark_repeats",1))
    wildcard_constraints:
        caller="(clairs_to|vardict|deepsomatic)"
    threads:
        config.get("bcftools_view",{}).get("threads",config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bcftools_view",{}).get("mem_mb",config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_view",{}).get("mem_per_cpu",config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_view",{}).get("partition",config["default_resources"]["partition"]),
        threads=config.get("bcftools_view",{}).get("threads",config["default_resources"]["threads"]),
        time=config.get("bcftools_view",{}).get("time",config["default_resources"]["time"]),
    container:
        config.get("bcftools_view",{}).get("container",config["default_container"])
    message:
        "{rule}: Use bcftools view to get subset or filter {input.vcf}"
    shell:
        """
        bcftools view {params.extra} -Oz -o {output.vcf} {input.vcf} 2> {log}
        """

rule filtering_bcftools_include_region:
    input:
        vcf="cnv_sv/{caller}/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.{caller}.vcf.gz",
        tbi="cnv_sv/{caller}/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.{caller}.vcf.gz.tbi",
    output:
        vcf=temp("cnv_sv/{caller}/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.{caller}.bcftools_view.vcf.gz"),
    params:
        extra=lambda wildcards: config.get("bcftools_view",{}) \
                                    .get("filter", {}) \
                                    .get(f"{wildcards.caller}",{}).get("extra", "")
    log:
        "cnv_sv/{caller}/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.{caller}.bcftools_view.vcf.log",
    benchmark:
        repeat("cnv_sv/{caller}/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.{caller}.bcftools_view.vcf.benchmark.tsv",
            config.get("bcftools_view",{}).get("benchmark_repeats",1))
    wildcard_constraints:
        caller="sniffles2"
    threads:
        config.get("bcftools_view",{}).get("threads",config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bcftools_view",{}).get("mem_mb",config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_view",{}).get("mem_per_cpu",config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_view",{}).get("partition",config["default_resources"]["partition"]),
        threads=config.get("bcftools_view",{}).get("threads",config["default_resources"]["threads"]),
        time=config.get("bcftools_view",{}).get("time",config["default_resources"]["time"]),
    container:
        config.get("bcftools_view",{}).get("container",config["default_container"])
    message:
        "{rule}: Use bcftools view to get subset or filter {input.vcf}"
    shell:
        """
        bcftools filter {params.extra} -Oz -o {output.vcf} {input.vcf} 2> {log}
        """


rule rename_vaf_to_af:
    input:
        vcf="{file}.vcf.gz",
        tbi="{file}.vcf.gz.tbi",
        txt="config/rename_vaf_to_af.txt"
    output:
        vcf=temp("{file}.rename_vaf.vcf"),
    wildcard_constraints:
        file="snv_indels/.*vep_annotated",
    log:
        "{file}.rename_vaf.vcf.gz.log",
    benchmark:
        repeat(
        "{file}.rename_vaf.vcf.gz.benchmark.tsv",
        config.get("vt_decompose",{}).get("benchmark_repeats",1),
        )
    resources:
        mem_mb=config.get("rename_vaf_to_af",{}).get("mem_mb",config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("rename_vaf_to_af",{}).get("mem_per_cpu",config["default_resources"]["mem_per_cpu"]),
        partition=config.get("rename_vaf_to_af",{}).get("partition",config["default_resources"]["partition"]),
        threads=config.get("rename_vaf_to_af",{}).get("threads",config["default_resources"]["threads"]),
        time=config.get("rename_vaf_to_af",{}).get("time",config["default_resources"]["time"]),
    container:
        config.get("rename_vaf_to_af",{}).get("container",config["default_container"])
    message:
        "{rule}: Use bcftools view to rename VAF to AF"
    shell:
        """
        bcftools annotate --force --rename-annots {input.txt} -o {output.vcf} {input.vcf} 2> {log}
        """