__author__ = "Camille Clouard"
__copyright__ = "Copyright 2024, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"

import os

# Include pipeline specific rules
include: "rules/common.smk"
include: "rules/basecalling.smk"
include: "rules/prefiltering.smk"
include: "rules/aligning.smk"
include: "rules/qc.smk"
include: "rules/varcalling.smk"
include: "rules/results_report.smk"

ruleorder: pycoqc > basecalling_dorado
ruleorder: fetch_filtered_reads > dorado_align
ruleorder: aligning_bam_index > mosdepth_overlap
ruleorder: mosdepth_overlap > qc_mosdepth_amplicons > varcall_clairs_to
ruleorder: filtlong > sequali
ruleorder: mosdepth_merge > qc_multiqc
ruleorder: mosdepth_overlap > varcall_clairs_to


# 'All' rule, must be specified before any other modules are
# included, since they also contain 'All' rule
rule all:
    input:
        "annotation/bcbio_variation_recall_ensemble/D22-06934_T.ensembled.vep_annotated.view.vcf"
        # unpack(compile_output_file_list),

# TODO: Draw DAG and rulegraph --> workflow/scripts
# snakemake -s workflow/Snakefile --configfile config/config.yaml --dag | dot -Tpdf > images/dag.pdf
# snakemake -s workflow/Snakefile --configfile config/config.yaml --rulegraph | dot -Tpdf > images/rulegraph.pdf

rule bcf_view_o_vcf:
    input:
        "annotation/bcbio_variation_recall_ensemble/D22-06934_T.ensembled.vep_annotated.vcf",
    output:
        "annotation/bcbio_variation_recall_ensemble/D22-06934_T.ensembled.vep_annotated.view.vcf.gz",
    params:
        extra="-Oz",
    wrapper:
        "v6.1.0/bio/bcftools/view"

# Include modules

module annotation:
    snakefile:
        github(
            "hydra-genetics/annotation",
            path="workflow/Snakefile",
            tag="v1.0.0",
        )
    config:
        config

# use rule * from annotation as annotation_*

use rule vep from annotation as annotation_vep with:
    input:
        cache=config.get("vep", {}).get("vep_cache", ""),
        fasta=config["ref_data"],
        tabix="snv_indels/clairs_to/{file}.vcf.gz.tbi",
        vcf="snv_indels/clairs_to/{file}.vcf.gz",
    output:
        vcf=temp("snv_indels/clairs_to/{file}.vep_annotated.vcf"),
    params:
        extra=config.get("vep", {}).get("extra", "--pick"),
        mode=config.get("vep", {}).get("mode", "--offline --cache --refseq "),
    log:
        "snv_indels/clairs_to/{file}.vep_annotated.vcf.log",
    benchmark:
        repeat("snv_indels/clairs_to/{file}.vep_annotated.vcf.benchmark.tsv",
            config.get("vep", {}).get("benchmark_repeats", 1))

# use rule bgzip_vcf from annotation as annotation_bgzip_vcf with:
#     input:
#         vcf=temp("snv_indels/clairs_to/{file}.vep_annotated.vcf"),
#     output:
#         gz=temp("snv_indels/clairs_to/{file}.vep_annotated.vcf.gz"),
#     log:
#         "snv_indels/clairs_to/{file}.vep_annotated.vcf.gz.log",
#     benchmark:
#         repeat("snv_indels/clairs_to/{file}.vep_annotated.vcf.gz.benchmark.tsv",
#             config.get("bgzip_vcf", {}).get("benchmark_repeats", 1))

# use rule tabix_vcf from annotation as annotation_tabix_vcf with:
#     input:
#         gz="snv_indels/clairs_to/{file}.vep_annotated.vcf.gz"
#     output:
#         tbi=temp("snv_indels/clairs_to/{file}.vep_annotated.vcf.gz.tbi")
#     log:
#         "snv_indels/clairs_to/{file}.vep_annotated.vcf.gz.tbi.log",
#     benchmark:
#         repeat(    "snv_indels/clairs_to/{file}.vep_annotated.vcf.gz.tbi.benchmark.tsv",
#             config.get("tabix_vcf",{}).get("benchmark_repeats",1))


module misc:
    snakefile:
        github(
            "hydra-genetics/misc",
            path="workflow/Snakefile",
            tag="v0.2.0",
        )
    config:
        config

# use rule * from misc as misc_*

use rule bgzip from misc as misc_bgzip with:
    input:
        gz="snv_indels/{file}.vcf",
    output:
        tbi=temp("snv_indels/{file}.vcf.gz"),
    params:
        extra=config.get("bgzip",{}).get("extra",""),
    log:
        "snv_indels/{file}.vcf.gz.log",
    benchmark:
        repeat(
            "snv_indels/{file}.vcf.gz.benchmark.tsv",
            config.get("bgzip",{}).get("benchmark_repeats",1),
        )

use rule tabix from misc as misc_tabix with:
    input:
        gz="snv_indels/{file}.vcf.gz",
    output:
        tbi=temp("snv_indels/{file}.vcf.gz.tbi"),
    params:
        extra=config.get("tabix",{}).get("extra",""),
    log:
        "snv_indels/{file}.vcf.gz.tbi.log",
    benchmark:
        repeat(
            "snv_indels/{file}.vcf.gz.tbi.benchmark.tsv",
            config.get("tabix",{}).get("benchmark_repeats",1),
        )


module qc:
    snakefile:
        github(
            "hydra-genetics/qc",
            path="workflow/Snakefile",
            tag=config["modules"]["qc"],
        )
    config:
        config

use rule * from qc as qc_*

use rule multiqc from qc as qc_multiqc with:
    input:
        files=lambda wildcards: set(
            [
                file.format(
                    sample=sample, target=target, type=type,
                )
                for file in config["multiqc"]["reports"][wildcards.report]["qc_files"]
                for sample in [wildcards.sample]
                for type in [wildcards.type]
                for target in config["amplicons"] + config["extra_regions"]
            ]
        ),
    output:
        html="results/qc/multiqc/{sample}_{type}_multiqc_{report}.html",
        data=directory("results/qc/multiqc/{sample}_{type}_multiqc_{report}_data")
    log:
        "results/qc/multiqc/{sample}_{type}_multiqc_{report}.html.log",
    benchmark:
        repeat("results/qc/multiqc/{sample}_{type}_multiqc_{report}.html.benchmark.tsv",
            config.get("multiqc",{}).get("benchmark_repeats",1))

# use rule samtools_stats from qc as qc_samtools_stats with:
#     input:
#         bam="{run_date}_{run_id}/data/align/reads.ont_adapt_trim.filtered.aligned.sorted.bam",
#     output:
#         stats="results/{run_date}_{run_id}/qc/samtools_stats/samtools-stats.txt",
#     params:
#         extra=f"--ref-seq {config['ref_data']} --target-regions {config['bed_files'] + '/amplicons.bed'}"
#     log:
#         "results/{run_date}_{run_id}/qc/samtools_stats/samtools-stats.txt.log"
#     benchmark:
#         repeat(
#             "results/{run_date}_{run_id}/qc/samtools_stats/samtools-stats.txt.benchmark.tsv",
#             config.get("samtools_stats",{}).get("benchmark_repeats",1),
#         )

use rule picard_collect_hs_metrics from qc as qc_picard_collect_hs_metrics with:
    input:
        bam="alignment/dorado_align/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam",
        bai="alignment/dorado_align/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam.bai",
        bait_intervals=config.get("reference", {}).get("design_intervals", ""),
        reference=config["ref_data"],
        target_intervals=config.get("reference", {}).get("design_intervals", ""),
    output:
        metrics = "results/qc/picard/{sample}_{type}_HsMetrics.txt",
    params:
        extra=config.get("picard_collect_hs_metrics", {}).get("extra", "COVERAGE_CAP=25000"),
    log:
        "results/qc/picard/{sample}_{type}_HsMetrics.txt.log",
    benchmark:
        repeat(
            "results/qc/picard/{sample}_{type}_HsMetrics.txt.benchmark.tsv",
            config.get("picard_collect_hs_metrics", {}).get("benchmark_repeats", 1),
        )


use rule mosdepth from qc as qc_mosdepth_amplicons with:
    input:
        bam="alignment/dorado_align/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam",
        bai="alignment/dorado_align/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam.bai",
    output:
        bed=temp("results/mosdepth/{sample}_{type}_amplicons.regions.bed.gz"),
        csi=temp("results/mosdepth/{sample}_{type}_amplicons.regions.bed.gz.csi"),
        glob=temp("results/mosdepth/{sample}_{type}_amplicons.mosdepth.global.dist.txt"),
        region=temp("results/mosdepth/{sample}_{type}_amplicons.mosdepth.region.dist.txt"),
        summary=temp("results/mosdepth/{sample}_{type}_amplicons.mosdepth.summary.txt"),
    params:
        by=os.path.join(config["bed_files"], config.get("mosdepth", {}).get("by", "")),
        extra=config.get("mosdepth", {}).get("extra", ""),
    log:
        "results/mosdepth/{sample}_{type}_amplicons.mosdepth.summary.txt.log",
    benchmark:
        repeat(
            "results/mosdepth/{sample}_{type}_amplicons.mosdepth.summary.txt.benchmark.tsv",
            config.get("mosdepth", {}).get("benchmark_repeats", 1),
        )
#     threads: config.get("mosdepth", {}).get("threads", config["default_resources"]["threads"])
#     resources:
#         mem_mb=config.get("mosdepth", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
#         mem_per_cpu=config.get("mosdepth", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
#         partition=config.get("mosdepth", {}).get("partition", config["default_resources"]["partition"]),
#         threads=config.get("mosdepth", {}).get("threads", config["default_resources"]["threads"]),
#         time=config.get("mosdepth", {}).get("time", config["default_resources"]["time"]),
#     container:
#         config.get("mosdepth", {}).get("container", config["default_container"])
