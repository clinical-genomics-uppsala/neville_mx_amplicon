__author__ = "Camille clouard"
__copyright__ = "Copyright 2026, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"


rule reference_rules_snv_indels_clairs_to_gvcf:
    input:
        bam="alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam",
        bai="alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam.bai",
        ref=config.get("ref_data"),
        bed=config.get("snv_indels_svs_clairs_to", {}).get("bed_file", os.path.join(config.get("bed_files"), "amplicons.bed")),
    output:
        snv=temp("references/clairs_to/{experiment}_{sample}_{type}/{experiment}_{sample}_{type}_snv.g.vcf.gz"),
        indel=temp("references/clairs_to/{experiment}_{sample}_{type}/{experiment}_{sample}_{type}_indel.g.vcf.gz"),
    params:
        platform=config.get("snv_indels_svs_clairs_to", {}).get("platform", ""),
        snv_min_af=config.get("snv_indels_svs_clairs_to", {}).get("snv_min_af", 0.05),
        indel_min_af=config.get("snv_indels_svs_clairs_to", {}).get("indel_min_af", 0.1),
        outdir=directory(lambda wildcards, output: os.path.dirname(output.snv)),
    resources:
        partition=config.get("snv_indels_svs_clairs_to", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("snv_indels_svs_clairs_to", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("snv_indels_svs_clairs_to", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("snv_indels_svs_clairs_to", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("snv_indels_svs_clairs_to", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("snv_indels_svs_clairs_to", {}).get("threads", config["default_resources"]["threads"])
    log:
        "references/clairs_to/{experiment}_{sample}_{type}_clairs_to.log",
    benchmark:
        repeat(
            "references/clairs_to/{experiment}_{sample}_{type}_clairs_to.benchmark.tsv",
            config.get("clairs_to", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("clairs_to", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Long-read somatic small variant calling in only tumor samples with ClairS-TO.
        """
    shell:
        "echo {params.outdir} && run_clairs_to --tumor_bam_fn {input.bam} --ref_fn {input.ref}"
        " --threads {resources.threads} --platform {params.platform}"
        " --output_dir {params.outdir} -s {wildcards.sample} --bed_fn {input.bed}"
        " --snv_min_af {params.snv_min_af} --indel_min_af {params.indel_min_af}"
        " --disable_verdict"
        " --print_ref_calls"
        " --snv_output_prefix {wildcards.experiment}_{wildcards.sample}_{wildcards.type}_snv"
        " --indel_output_prefix {wildcards.experiment}_{wildcards.sample}_{wildcards.type}_indel"
        " 2> {log}"

# rule reference_rules_create_background_file_clairs_to: