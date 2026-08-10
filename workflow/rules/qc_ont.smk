__author__ = "Camille Clouard"
__copyright__ = "Copyright 2024, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"

from snakemake.logging import logger

logger.info(f"\n{workflow.snakefile} is being parsed")


rule qc_ont_pycoqc:
    input:
        seq_run_dir=os.path.join(config.get("runfolder"), config.get("batchid"), config.get("runid")),
    output:
        html=temp("results/pycoqc/{experiment}_{sample}_{type}_report_sequencing_summary.html"),
        json=temp("results/pycoqc/{experiment}_{sample}_{type}_report_sequencing_summary.json"),
        txt=temp("results/pycoqc/{experiment}_{sample}_{type}_report_sequencing_summary.txt"),
    resources:
        partition=config.get("qc_ont_pycoqc", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("qc_ont_pycoqc", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("qc_ont_pycoqc", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("qc_ont_pycoqc", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("qc_ont_pycoqc", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("qc_ont_pycoqc", {}).get("threads", config["default_resources"]["threads"])
    log:
        "results/pycoqc/{experiment}_{sample}_{type}_report_sequencing_summary.log",
    benchmark:
        repeat(
            "results/pycoqc/{experiment}_{sample}_{type}_report_sequencing_summary.benchmark.tsv",
            config.get("qc_ont_pycoqc", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("qc_ont_pycoqc", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Report graphically for the sequencing run.
        """
    shell:
        """
        # summary=$( ls {input.seq_run_dir}/sequencing_summary*.txt )  # fails in some shells
        summary=$(find {input.seq_run_dir} -maxdepth 1 -name "sequencing_summary*.txt" | head -n 1)
        cp $summary {output.txt}
        pycoQC -f {output.txt} --html_outfile {output.html} --json_outfile {output.json} 2> {log}
        """


rule qc_ont_mosdepth_overlap:
    input:
        bam="alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam",
        bamidx="alignment/dorado_align/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam.bai",
        amplibed=os.path.join(config.get("bed_files"), "{target}.bed"),
    output:
        bed=temp("results/mosdepth/{experiment}_{sample}_{type}_{target}.regions.bed.gz"),
        csi=temp("results/mosdepth/{experiment}_{sample}_{type}_{target}.regions.bed.gz.csi"),
        glob=temp("results/mosdepth/{experiment}_{sample}_{type}_{target}.mosdepth.global.dist.txt"),
        region=temp("results/mosdepth/{experiment}_{sample}_{type}_{target}.mosdepth.region.dist.txt"),
        summary=temp("results/mosdepth/{experiment}_{sample}_{type}_{target}.mosdepth.summary.txt"),
    params:
        prefix_out=lambda wildcards, output: os.path.dirname(output.summary),
        threads=20,
    resources:
        partition=config.get("qc_ont_mosdepth_overlap", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("qc_ont_mosdepth_overlap", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("qc_ont_mosdepth_overlap", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("qc_ont_mosdepth_overlap", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("qc_ont_mosdepth_overlap", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("qc_ont_mosdepth_overlap", {}).get("threads", config["default_resources"]["threads"])
    log:
        "results/mosdepth/{experiment}_{sample}_{type}_{target}.mosdepth.log",
    benchmark:
        repeat(
            "results/mosdepth/{experiment}_{sample}_{type}_{target}.mosdepth.benchmark.tsv",
            config.get("qc_ont_mosdepth_overlap", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("qc_ont_mosdepth_overlap", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Compute coverage with mosdepth for each amplicon.
        """
    shell:
        """
        chrom=$( cat {input.amplibed} | cut -d$'\t' -f1 )
        mosdepth -t {params.threads} -c $chrom -b {input.amplibed}  {params.prefix_out}/{wildcards.experiment}_{wildcards.sample}_{wildcards.type}_{wildcards.target} {input.bam} > {log}
        """


rule qc_ont_mosdepth_merge:
    input:
        expand(
            "results/mosdepth/{{experiment}}_{{sample}}_{{type}}_{target}.mosdepth.summary.txt",
            target=config.get("amplicons") + config.get("extra_regions"),
        ),
    output:
        csv=temp("results/mosdepth/{experiment}_{sample}_{type}_coverage_per_amplicon.csv"),
    resources:
        partition=config.get("default_resources").get("partition"),
        time=config.get("default_resources").get("time"),
        threads=config.get("default_resources").get("threads"),
        mem_mb=config.get("default_resources").get("mem_mb"),
        mem_per_cpu=config.get("default_resources").get("mem_per_cpu"),
    threads: config.get("default_resources").get("threads")
    log:
        "results/mosdepth/{experiment}_{sample}_{type}_coverage_per_amplicon.log",
    benchmark:
        repeat(
            "results/mosdepth/{experiment}_{sample}_{type}_coverage_per_amplicon.benchmark.tsv",
            config.get("qc_ont_mosdepth_merge", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("qc_ont_mosdepth_merge", {}).get("container", config["default_container"])
    message:
        "{rule}: Create merged report for mosdepth"
    script:
        "../scripts/mosdepth_merge.py"


rule qc_ont_mosdepth_overlap_timestep:
    input:
        bamdir=os.path.join(config["runfolder"], "{sample}", config["runid"], "bam_pass"),
        amplibed=[f"{config.get('bed_files')}/{target}.bed" for target in config.get("amplicons") + config.get("extra_regions")],
    output:
        outdir=temp(directory("results/mosdepth/timestep/{experiment}_{sample}")),
    resources:
        partition=config.get("qc_ont_mosdepth_overlap_timestep", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("qc_ont_mosdepth_overlap_timestep", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("qc_ont_mosdepth_overlap_timestep", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("qc_ont_mosdepth_overlap_timestep", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("qc_ont_mosdepth_overlap_timestep", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
    threads: config.get("qc_ont_mosdepth_overlap_timestep", {}).get("threads", config["default_resources"]["threads"])
    log:
        "results/mosdepth/timestep/{experiment}_{sample}.log",
    benchmark:
        repeat(
            "results/mosdepth/timestep/{experiment}_{sample}.tsv",
            config.get("qc_ont_mosdepth_overlap_timestep", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("qc_ont_mosdepth_overlap_timestep", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Compute coverage with mosdepth for each amplicon.
        """
    script:
        "../scripts/process_timestep_data.py"


rule qc_ont_mosdepth_merge_timestep:
    input:
        indir="results/mosdepth/timestep/{experiment}_{sample}",
    output:
        outdir=temp(directory("results/mosdepth/timestep_coverage/{experiment}_{sample}")),
    resources:
        partition=config.get("default_resources").get("partition"),
        time=config.get("default_resources").get("time"),
        threads=config.get("default_resources").get("threads"),
        mem_mb=config.get("default_resources").get("mem_mb"),
        mem_per_cpu=config.get("default_resources").get("mem_per_cpu"),
    threads: config.get("default_resources").get("threads")
    container:
        config.get("qc_ont_mosdepth_merge_timestep", {}).get("container", config["default_container"])
    log:
        "results/mosdepth/timestep_coverage/{experiment}_{sample}/timestep_coverage.log",
    benchmark:
        repeat(
            "results/mosdepth/timestep_coverage/{experiment}_{sample}/timestep_coverage.benchmark.tsv",
            config.get("qc_ont_mosdepth_merge_timestep", {}).get("benchmark_repeats", 1),
        )
    message:
        "{rule}: Create merged report for mosdepth"
    script:
        "../scripts/mosdepth_merge_timestep.py"


rule qc_ont_plot_yield_timestep:
    input:
        indir="results/mosdepth/timestep_coverage/{experiment}_{sample}",
    output:
        csv=temp("results/mosdepth/timestep_coverage_images/{experiment}_{sample}_{type}_cumsum_coverage_per_amplicon.csv"),
        png=temp("results/mosdepth/timestep_coverage_images/{experiment}_{sample}_{type}_cumsum_coverage_per_amplicon.png"),
    params:
        timestep=config.get("timestep_minknow", 10)
    resources:
        partition=config["default_resources"]["partition"],
        time=config["default_resources"]["time"],
        threads=config["default_resources"]["threads"],
        mem_mb=config["default_resources"]["mem_mb"],
        mem_per_cpu=config["default_resources"]["mem_per_cpu"],
    threads: config["default_resources"]["threads"]
    container:
        config.get("qc_ont_plot_yield_timestep", {}).get("container", config["default_container"])
    log:
        "results/mosdepth/timestep_coverage_images/{experiment}_{sample}_{type}_cumsum_coverage_per_amplicon.log",
    benchmark:
        repeat(
            "results/mosdepth/timestep_coverage_images/{experiment}_{sample}_{type}_cumsum_coverage_per_amplicon.benchmark.tsv",
            config.get("qc_ont_plot_yield_timestep", {}).get("benchmark_repeats", 1),
        )
    message:
        """
        {rule}: Plot sequencing output per amplicon over time.
        """
    script:
        "../scripts/seq_yield_timestep.py"


rule qc_ont_sequali:
    input:
        fastgz1="prealignment/filtlong/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.fastq.gz",
        fastgz2="prealignment/filtlong/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.out.fastq.gz",
    output:
        html1=temp("results/sequali/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.fastq.gz.html"),
        json1=temp("results/sequali/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.fastq.gz.json"),
        html2=temp("results/sequali/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.out.fastq.gz.html"),
        json2=temp("results/sequali/{experiment}_{sample}_{type}_reads.ont_adapt_trim.filtered.out.fastq.gz.json"),
    resources:
        partition=config.get("qc_ont_sequali", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("qc_ont_sequali", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("qc_ont_sequali", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("qc_ont_sequali", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("qc_ont_sequali", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("qc_ont_sequali", {}).get("threads", config["default_resources"]["threads"])
    log:
        "results/sequali/{experiment}_{sample}_{type}_sequali.log",
    benchmark:
        repeat(
            "results/sequali/{experiment}_{sample}_{type}_sequali.benchmark.tsv",
            config.get("qc_ont_sequali", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("qc_ont_sequali", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Parse the (out)filtered reads and create quality report with sequali.
        """
    shell:
        """
        sequali --html {output.html1} --json {output.json1} {input.fastgz1} 2> {log}
        sequali --html {output.html2} --json {output.json2} {input.fastgz2} 2>> {log}
        """


rule qc_ont_yield_per_pool:
    input:
        csv="results/mosdepth/{experiment}_{sample}_{type}_coverage_per_amplicon.csv",
    output:
        csv=temp("results/mosdepth/{experiment}_{sample}_{type}_yield_pool_{pooln}.csv"),
    resources:
        partition=config.get("default_resources").get("partition"),
        time=config.get("default_resources").get("time"),
        threads=config.get("default_resources").get("threads"),
        mem_mb=config.get("default_resources").get("mem_mb"),
        mem_per_cpu=config.get("default_resources").get("mem_per_cpu"),
    threads: config.get("default_resources").get("threads")
    container:
        config.get("qc_ont_yield_per_pool", {}).get("container", config["default_container"])
    log:
        "results/mosdepth/{experiment}_{sample}_{type}_yield_pool_{pooln}.log",
    benchmark:
        repeat(
            "results/mosdepth/{experiment}_{sample}_{type}_yield_pool_{pooln}.benchmark.tsv",
            config.get("qc_ont_yield_per_pool", {}).get("benchmark_repeats", 1),
        )
    message:
        "{rule}: Calculate number of reads per pool"
    script:
        "../scripts/yield_per_pool.py"


# Picard HsMetrics requires interval file
rule qc_ont_picard_bed_to_interval_list:
    input:
        bed=os.path.join(config["bed_files"], "amplicons.bed"),
        dict=config["reference"]["sequence_dict"],
    output:
        config["reference"]["design_intervals"],
    log:
        "results/qc/picard/BedToIntervalList.log",
    params:
        extra="--SORT true",  # sort output interval list before writing
    resources:
        partition=config.get("default_resources").get("partition"),
        time=config.get("default_resources").get("time"),
        threads=config.get("default_resources").get("threads"),
        mem_mb=config.get("default_resources").get("mem_mb"),
        mem_per_cpu=config.get("default_resources").get("mem_per_cpu"),
    threads: config.get("default_resources").get("threads")
    benchmark:
        repeat(
            "results/qc/picard/BedToIntervalList.benchmark.tsv",
            config.get("qc_ont_picard_bed_to_interval_list", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("qc_ont_picard_bed_to_interval_list", {}).get("container", config["default_container"])
    wrapper:
        "v5.0.1/bio/picard/bedtointervallist"


rule qc_ont_mosdepth_exons_multiqc:
    input:
        bed="results/mosdepth_bed_per_exon/{experiment}_{sample}_{type}.regions.bed.gz",
    output:
        tsv="results/qc/custom_content_mqc/{experiment}_{sample}_{type}.exon_coverage_mqc.tsv",
    params:
        extra=config.get("qc_ont_mosdepth_exons_multiqc", {}).get("extra", ""),
    log:
        "results/qc/custom_content_mqc/{experiment}_{sample}_{type}.exon_coverage_mqc.output.log",
    benchmark:
        repeat(
            "results/qc/custom_content_mqc/{experiment}_{sample}_{type}.exon_coverage_mqc.output.benchmark.tsv",
            config.get("qc_ont_mosdepth_exons_multiqc", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("qc_ont_mosdepth_exons_multiqc", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("qc_ont_mosdepth_exons_multiqc", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("qc_ont_mosdepth_exons_multiqc", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("qc_ont_mosdepth_exons_multiqc", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("qc_ont_mosdepth_exons_multiqc", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("qc_ont_mosdepth_exons_multiqc", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("qc_ont_mosdepth_exons_multiqc", {}).get("container", config["default_container"])
    message:
        "{rule}: Write coverage values in exons to TSV file"
    script:
        "../scripts/mosdepth_exons_multiqc.py"
