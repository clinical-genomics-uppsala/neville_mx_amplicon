__author__ = "Camille clouard"
__copyright__ = "Copyright 2026, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"


rule reference_rules_snv_indels_clairs_to_gvcf:
    input:
        bam="alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam",
        bai="alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam.bai",
        ref=config.get("ref_data"),
        bed=config.get("reference_rules_snv_indels_clairs_to_gvcf", {}).get("bed_file", os.path.join(config.get("bed_files"), "amplicons.bed")),
        vcf="references/deepsomatic/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.deepsomatic.g.vcf.gz",
    output:
        snv="references/clairs_to/{experiment}_{sample}_{type}/{experiment}_{sample}_{type}_snv.g.vcf.gz",
        # indel="references/clairs_to/{experiment}_{sample}_{type}/{experiment}_{sample}_{type}_indel.g.vcf.gz",  # no indel file is written for some reason: is it just for that sample or in general when using the -H/-G option?
    params:
        platform=config.get("reference_rules_snv_indels_clairs_to_gvcf", {}).get("platform", ""),
        snv_min_af=config.get("reference_rules_snv_indels_clairs_to_gvcf", {}).get("snv_min_af", 0.05),
        indel_min_af=config.get("reference_rules_snv_indels_clairs_to_gvcf", {}).get("indel_min_af", 0.1),
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
            config.get("reference_rules_snv_indels_clairs_to_gvcf", {}).get("benchmark_repeats", 1),
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
        # " -G {input.vcf}"
        " --print_ref_calls"
        " --snv_output_prefix {wildcards.experiment}_{wildcards.sample}_{wildcards.type}_snv.g"
        " --indel_output_prefix {wildcards.experiment}_{wildcards.sample}_{wildcards.type}_indel.g"
        " 2> {log}"


rule reference_rules_snv_indels_deepsomatic_gvcf:
    input:
        bam="alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam",
        bai="alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam.bai",
        ref=config.get("ref_data"),
        bed=config.get("reference_rules_snv_indels_deepsomatic_gvcf", {}).get("bed_file", os.path.join(config.get("bed_files"), "amplicons.bed")),
    output:
        tmpdir=directory("references/deepsomatic/{experiment}_{sample}_{type}_tmp"),
        vcf=temp("references/deepsomatic/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.deepsomatic.vcf.gz"),
        gvcf=temp("references/deepsomatic/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.deepsomatic.g.vcf.gz"),
    params:
        sample=lambda wildcards: f"{wildcards.sample}",
        model="ONT_TUMOR_ONLY",
        extra=config.get("reference_rules_snv_indels_deepsomatic_gvcf", {}).get("extra", ""),
        filter=config.get("reference_rules_snv_indels_deepsomatic_gvcf", {}).get("filter", ""),
    log:
        "references/deepsomatic/{experiment}_{sample}_{type}_deepsomatic.log",
    benchmark:
        repeat(
            "references/deepsomatic/{experiment}_{sample}_{type}_deepsomatic.benchmark.tsv",
            config.get("reference_rules_snv_indels_deepsomatic_gvcf", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("reference_rules_snv_indels_deepsomatic_gvcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        partition=config.get("snv_indels_svs_deepsomatic", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("snv_indels_svs_deepsomatic", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("snv_indels_svs_deepsomatic", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("snv_indels_svs_deepsomatic", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("snv_indels_svs_deepsomatic", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    container:
        config.get("reference_rules_snv_indels_deepsomatic_gvcf", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Long-read somatic small variant calling in tumor samples only with DeepSomatic.
        """
    shell:
        """
        run_deepsomatic --model_type={params.model} --ref={input.ref} --reads_tumor={input.bam} --output_vcf={output.vcf} --output_gvcf={output.gvcf} --sample_name_tumor={params.sample} --num_shards={resources.threads} --logging_dir={log} --intermediate_results_dir {output.tmpdir} --regions={input.bed} {params.extra} {params.filter}
        """


rule create_background_file_longread_clairs_to:
    input:
        gvcfs=get_gvcfs_clairs_to(),
    output:
        background_file=temp("references/create_background_file/background_panel_clairs_to.tsv"),
    params:
        min_dp=config.get("create_background_file", {}).get("min_dp", 300),
        max_af=config.get("create_background_file", {}).get("max_af", 0.03),
    log:
        "references/create_background_file/background_panel_clairs_to.tsv.log",
    benchmark:
        repeat(
            "references/create_background_file/background_panel_clairs_to.tsv.benchmark.tsv",
            config.get("create_background_file", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("create_background_file", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("create_background_file", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_background_file", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_background_file", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_background_file", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_background_file", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("create_background_file", {}).get("container", config["default_container"])
    message:
        "{rule}: create background PoN"
    script:
        "../scripts/create_background_file_longread.py"


rule create_background_file_longread_deepsomatic:
    input:
        gvcfs=get_gvcfs_deepsomatic(),
    output:
        background_file=temp("references/create_background_file/background_panel_deepsomatic.tsv"),
    params:
        min_dp=config.get("create_background_file", {}).get("min_dp", 300),
        max_af=config.get("create_background_file", {}).get("max_af", 0.03),
    log:
        "references/create_background_file/background_panel_deepsomatic.tsv.log",
    benchmark:
        repeat(
            "references/create_background_file/background_panel_deepsomatic.tsv.benchmark.tsv",
            config.get("create_background_file", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("create_background_file", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("create_background_file", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_background_file", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_background_file", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_background_file", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_background_file", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("create_background_file", {}).get("container", config["default_container"])
    message:
        "{rule}: create background PoN"
    script:
        "../scripts/create_background_file_longread.py"


rule create_background_file_longread_vardict:
    input:
        gvcfs=get_gvcfs_vardict(),
    output:
        background_file=temp("references/create_background_file/background_panel_vardict.tsv"),
    params:
        min_dp=config.get("create_background_file", {}).get("min_dp", 300),
        max_af=config.get("create_background_file", {}).get("max_af", 0.03),
    log:
        "references/create_background_file/background_panel_vardict.tsv.log",
    benchmark:
        repeat(
            "references/create_background_file/background_panel_vardict.tsv.benchmark.tsv",
            config.get("create_background_file", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("create_background_file", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("create_background_file", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_background_file", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_background_file", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_background_file", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_background_file", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("create_background_file", {}).get("container", config["default_container"])
    message:
        "{rule}: create background PoN"
    script:
        "../scripts/create_background_file_longread.py"


rule create_background_file_longread_vardict_unfiltered:
    input:
        gvcfs=get_gvcfs_vardict_unfiltered(),
    output:
        background_file=temp("references/create_background_file/background_panel_vardict_unfiltered.tsv"),
    params:
        min_dp=config.get("create_background_file", {}).get("min_dp", 300),
        max_af=config.get("create_background_file", {}).get("max_af", 0.03),
    log:
        "references/create_background_file/background_panel_vardict_unfiltered.tsv.log",
    benchmark:
        repeat(
            "references/create_background_file/background_panel_vardict_unfiltered.tsv.benchmark.tsv",
            config.get("create_background_file", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("create_background_file", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("create_background_file", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_background_file", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_background_file", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_background_file", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_background_file", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("create_background_file", {}).get("container", config["default_container"])
    message:
        "{rule}: create background PoN"
    script:
        "../scripts/create_background_file_longread.py"


rule create_artifact_file_longread:
    input:
        vcfs=get_vcfs(),
    output:
        artifact_panel=temp("references/create_artifact_file/artifact_panel.tsv"),
    params:
        callers=config.get("create_artifact_file", {}).get("callers", ["vardict", "deepsomatic", "clairs_to"]),
    log:
        "references/create_artifact_file/artifact_panel.tsv.log",
    benchmark:
        repeat(
            "references/create_artifact_file/artifact_panel.tsv.benchmark.tsv",
            config.get("create_artifact_file", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("create_artifact_file", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("create_artifact_file", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_artifact_file", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_artifact_file", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_artifact_file", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_artifact_file", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("create_artifact_file", {}).get("container", config["default_container"])
    message:
        "{rule}: create artifact PoN"
    script:
        "../scripts/create_artifact_file_longread.py"


rule reference_rules_align_unfiltered:
    input:
        fastqgz="basecalling/dorado_duplex/{experiment}_{sample}_{type}_reads.ont_adapt_trim.fastq.gz",
        ref_data=config.get("ref_data"),
    output:
        bam=temp("references/align_unfiltered_bam/{experiment}_{sample}_{type}_reads.ont_adapt_trim.unfiltered.aligned.bam"),
    params:
        runid=config.get("runid", ""),
        extra=config.get("alignment_ont_dorado_align", {}).get("extra", ""),
    resources:
        partition=config.get("alignment_ont_dorado_align", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("alignment_ont_dorado_align", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("alignment_ont_dorado_align", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("alignment_ont_dorado_align", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("alignment_ont_dorado_align", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("alignment_ont_dorado_align", {}).get("threads", config["default_resources"]["threads"])
    container:
        config.get("alignment_ont_dorado_align", {}).get("container", config["default_container"])
    log:
        "references/align_unfiltered_bam/{experiment}_{sample}_{type}_reads.ont_adapt_trim.unfiltered.aligned.log",
    benchmark:
        repeat(
            "references/align_unfiltered_bam/{experiment}_{sample}_{type}_reads.ont_adapt_trim.unfiltered.aligned.benchmark.tsv",
            config.get("alignment_ont_dorado_align", {}).get("benchmark_repeats", 1),
        )
    message:
        "{rule}: Align reads with dorado and minimap2"
    shell:
        """
        echo "Dorado executed from $( which dorado )"

        echo "Executing dorado aligning of {input.fastqgz} with reference genome '{input.ref_data}'"

        dorado aligner {params.extra} {input.ref_data} {input.fastqgz} > {output.bam} 2> {log}
        """


rule reference_rules_bam_sort_unfiltered:
    input:
        "references/align_unfiltered_bam/{experiment}_{sample}_{type}_reads.ont_adapt_trim.unfiltered.aligned.bam",
    output:
        temp("references/align_unfiltered_bam/{experiment}_{sample}_{type}_reads.ont_adapt_trim.unfiltered.aligned.sorted.bam"),
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
        "references/align_unfiltered_bam/{experiment}_{sample}_{type}_reads.ont_adapt_trim.unfiltered.aligned.sorted.log",
    benchmark:
        repeat(
            "references/align_unfiltered_bam/{experiment}_{sample}_{type}_reads.ont_adapt_trim.unfiltered.aligned.sorted.benchmark.tsv",
            config.get("alignment_ont_bam_sort", {}).get("benchmark_repeats", 1),
        )
    message:
        "{rule}: Sort aligned reads with samtools"
    wrapper:
        "0.2.0/bio/samtools/sort"


rule reference_rules_bam_index_unfiltered:
    input:
        "references/align_unfiltered_bam/{experiment}_{sample}_{type}_reads.ont_adapt_trim.unfiltered.aligned.sorted.bam",
    output:
        temp("references/align_unfiltered_bam/{experiment}_{sample}_{type}_reads.ont_adapt_trim.unfiltered.aligned.sorted.bam.bai"),
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
        "references/align_unfiltered_bam/{experiment}_{sample}_{type}_reads.ont_adapt_trim.unfiltered.aligned.sorted.bam.bai.log",
    benchmark:
        repeat(
            "references/align_unfiltered_bam/{experiment}_{sample}_{type}_reads.ont_adapt_trim.unfiltered.aligned.sorted.bam.bai.benchmark.tsv",
            config.get("alignment_ont_bam_index", {}).get("benchmark_repeats", 1),
        )
    message:
        "{rule}: Index the aligned and sorted reads with samtools"
    wrapper:
        "0.2.0/bio/samtools/index"