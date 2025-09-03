__author__ = "Camille Clouard"
__copyright__ = "Copyright 2025, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"

import os.path

from hydra_genetics.utils.software_versions import get_pipeline_version


# rule results_report_bedtools_intersect:
#     input:
#         left="qc/mosdepth_bed_coding/{sample}_{type}.per-base.bed.gz",
#         coverage_csi="qc/mosdepth_bed_coding/{sample}_{type}.per-base.bed.gz.csi",
#         right=config["reference"]["exon_bed"],
#     output:
#         results=temp("qc/mosdepth_bed_coding/{sample}_{type}.mosdepth.per-base.exon_bed.txt"),
#     params:
#         extra=config.get("results_report_bedtools_intersect", {}).get("extra", ""),
#     log:
#         "qc/mosdepth_bed_coding/{sample}_{type}.mosdepth.per-base.exon_bed.log",
#     benchmark:
#         repeat(
#             "qc/mosdepth_bed_coding/{sample}_{type}.mosdepth.per-base.exon_bed.benchmark.tsv",
#             config.get("results_report_bedtools_intersect", {}).get("benchmark_repeats", 1),
#         )
#     threads: config.get("results_report_bedtools_intersect", {}).get("threads", config["default_resources"]["threads"])
#     resources:
#         mem_mb=config.get("results_report_bedtools_intersect", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
#         mem_per_cpu=config.get("results_report_bedtools_intersect", {}).get(
#             "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
#         ),
#         partition=config.get("results_report_bedtools_intersect", {}).get("partition", config["default_resources"]["partition"]),
#         threads=config.get("results_report_bedtools_intersect", {}).get("threads", config["default_resources"]["threads"]),
#         time=config.get("results_report_bedtools_intersect", {}).get("time", config["default_resources"]["time"]),
#     container:
#         config.get("results_report_bedtools_intersect", {}).get("container", config["default_container"])
#     message:
#         "{rule}: export low cov regions from {input.left} based on {input.right}"
#     wrapper:
#         "v1.32.0/bio/bedtools/intersect"


rule results_report_xlsx:
    input:
        vcf_snv_indels="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.rename_vaf.filter.somatic_soft.vcf.gz",
        vcf_snv_indels_tbi="snv_indels/bcbio_variation_recall_ensemble/{sample}_{type}.ensembled.vep_annotated.rename_vaf.filter.somatic_soft.vcf.gz.tbi",
        # TODO: replace the following with tab for Sniffles2
        vcf_cnv_svs="cnv_sv/sniffles2/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.sniffles2.bcftools_view.vep_annotated.vcf.gz",
        vcf_cnv_svs_tbi="cnv_sv/sniffles2/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.sniffles2.bcftools_view.vep_annotated.vcf.gz.tbi",
        bedfile=os.path.join(config["bed_files"], "amplicons.bed"),
        deepsomatic_bed=config["deepsomatic"]["bed_file"],
        vardict_bed=config["vardict"]["bed_file"],
        # extra_bedfile=[config[caller]["bed_file"] for caller in ["deepsomatic", "vardict", "sniffles2"]
        mosdepth_summary="results/mosdepth/{sample}_{type}_amplicons.mosdepth.summary.txt",
        mosdepth_exons="results/mosdepth_bed_per_exon/{sample}_{type}.regions.bed.gz",
        mosdepth_regions="results/mosdepth/{sample}_{type}_amplicons.regions.bed.gz",
        csv_counts="results/mosdepth/{sample}_{type}_coverage_per_amplicon.csv",
        pool_counts=expand("results/mosdepth/{{sample}}_{{type}}_yield_pool_{pooln}.csv",
            pooln=config["pools"].keys()
        ),
        yield_plot="results/mosdepth/timestep_coverage_images/{sample}_{type}_cumsum_coverage_per_amplicon.png"
        # mosdepth_thresholds="qc/mosdepth_bed_coding/{sample}_{type}.thresholds.bed.gz",
        # wanted_transcripts=config["results_report_xlsx"]["wanted_transcripts"],
    output:
        xlsx="reports/xlsx/{sample}_{type}.xlsx",
    params:
        sample=lambda wildcards: wildcards.sample,
        sample_type=lambda wildcards: wildcards.type,
        sequenceid=config["runfolder"],
        # poppy_version=config["poppy_version"],
        # uppsala_version=get_pipeline_version(workflow, pipeline_name="poppy_uppsala"),
        ref=config["reference"]["fasta"],
        amplicons=config["amplicons"],
        extra_regions=config["extra_regions"],
        # artifact=config["reference"]["artifacts"],
        # background=config["reference"]["background"],
        # artifact_indel=config["reference"]["artifacts_indel"],
        # thresholds=config["mosdepth_bed"]["thresholds"],
        filter_somatic_soft=config["filter_vcf"]["somatic_soft"],
        filter_somatic_hard=config["filter_vcf"]["somatic_hard"],
        # filter_indel=config["filter_vcf"]["indel"],
    log:
        "reports/xlsx/{sample}_{type}.xlsx.log",
    benchmark:
        repeat("reports/xlsx/{sample}_{type}.xlsx.benchmark.tsv", config.get("results_report", {}).get("benchmark_repeats", 1))
    threads: config.get("results_report", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("results_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("results_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("results_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("results_report", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("results_report", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("results_report", {}).get("container", config["default_container"])
    message:
        "{rule}: summarize results into {output.xlsx}."
    script:
        "../scripts/results_report_xlsx.py"

rule copy_bed:
    input:
        [bed_file for  bed_file in [config[caller]["bed_file"] for caller in ["deepsomatic", "vardict"]]]+ [os.path.join(config["bed_files"], "amplicons.bed")]
    output:
        outdir=directory("bedfiles")
    log:
        "logs/copy_bed.log"
    resources:
        partition=config.get("copy_bed",{}).get("partition",config["default_resources"]["partition"]),
        time=config.get("copy_bed",{}).get("time",config["default_resources"]["time"]),
        threads=config.get("copy_bed",{}).get("threads",config["default_resources"]["threads"]),
        mem_mb=config.get("copy_bed",{}).get("mem_mb",config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("copy_bed",{}).get("mem_per_cpu",config["default_resources"]["mem_per_cpu"]),
    container:
        config.get("copy_bed",{}).get("container",config["default_container"])
    message:
        """
        {rule}: Copy the different BED files used in the analysis.
        """
    shell:
        """
        mkdir -p bedfiles
        for bed in {input}
        do
            echo $( basename $bed )
            cp $bed {output.outdir}/$( basename $bed ) 2>> {log}
        done
        """