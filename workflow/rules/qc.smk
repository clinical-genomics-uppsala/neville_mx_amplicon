__author__ = "Camille Clouard"
__copyright__ = "Copyright 2024, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"


def get_bam_pass_sample():
    expr = lambda wildcards: expand("results/mosdepth/timestep/{{fname}}_{{nbatch}}/{target}.mosdepth.summary.txt",
            fname=read_bam_pass_names(config["runfolder"], f"{wildcards.sample}", config["runid"], "bam_pass")[0],
            nbatch=read_bam_pass_names(config["runfolder"], f"{wildcards.sample}", config["runid"], "bam_pass")[1],
            target=config.get("amplicons") + config.get("extra_regions"),
        )
    return expr


rule pycoqc:
    input:
        seq_run_dir=os.path.join(config.get("runfolder"), config.get("batchid"), config.get("runid"))
    output:
        html=temp("results/pycoqc/{sample}_{type}_report_sequencing_summary.html"),
        json=temp("results/pycoqc/{sample}_{type}_report_sequencing_summary.json"),
        txt=temp("results/pycoqc/{sample}_{type}_report_sequencing_summary.txt"),
    resources:
        partition=config.get("pycoqc", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("pycoqc", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("pycoqc", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("pycoqc", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pycoqc", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("pycoqc", {}).get("threads", config["default_resources"]["threads"])
    log:
        "results/pycoqc/{sample}_{type}_report_sequencing_summary.log",
    benchmark:
        repeat("results/pycoqc/{sample}_{type}_report_sequencing_summary.benchmark.tsv",
            config.get("pycoqc", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("pycoqc", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Report graphically for the sequencing run.
        """
    shell:
        """
        summary=$( ls {input.seq_run_dir}/sequencing_summary*.txt ) 
        cp $summary {output.txt}
        pycoQC -f {output.txt} --html_outfile {output.html} --json_outfile {output.json} 2> {log}
        """


rule mosdepth_overlap:
    input:
        bam="alignment/dorado_align/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam",
        bamidx="alignment/dorado_align/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam.bai",
        amplibed=os.path.join(config.get("bed_files"), "{target}.bed"),
    output:
        bed=temp("results/mosdepth/{sample}_{type}_{target}.regions.bed.gz"),
        csi=temp("results/mosdepth/{sample}_{type}_{target}.regions.bed.gz.csi"),
        glob=temp("results/mosdepth/{sample}_{type}_{target}.mosdepth.global.dist.txt"),
        region=temp("results/mosdepth/{sample}_{type}_{target}.mosdepth.region.dist.txt"),
        summary=temp("results/mosdepth/{sample}_{type}_{target}.mosdepth.summary.txt"),
    params:
        prefix_out=lambda wildcards, output: os.path.dirname(output.summary),
        threads=20,
    resources:
        partition=config.get("mosdepth", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("mosdepth", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("mosdepth", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("mosdepth", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mosdepth", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("mosdepth", {}).get("threads", config["default_resources"]["threads"])
    log:
        "results/mosdepth/{sample}_{type}_{target}.mosdepth.log",
    benchmark:
        repeat("results/mosdepth/{sample}_{type}_{target}.mosdepth.benchmark.tsv",
            config.get("mosdepth", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("mosdepth", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Compute coverage with mosdepth for each amplicon.
        """
    shell:
        """
        chrom=$( cat {input.amplibed} | cut -d$'\t' -f1 )
        mosdepth -t {params.threads} -c $chrom -b {input.amplibed}  {params.prefix_out}/{wildcards.sample}_{wildcards.type}_{wildcards.target} {input.bam} 2> {log}
        """


rule mosdepth_merge:
    input:
        expand("results/mosdepth/{{sample}}_{{type}}_{target}.mosdepth.summary.txt",
            target=config.get("amplicons") + config.get("extra_regions"),
        ),
    output:
        csv=temp("results/mosdepth/{sample}_{type}_coverage_per_amplicon.csv"),
    resources:
        partition=config.get("default_resources").get("partition"),
        time=config.get("default_resources").get("time"),
        threads=config.get("default_resources").get("threads"),
        mem_mb=config.get("default_resources").get("mem_mb"),
        mem_per_cpu=config.get("default_resources").get("mem_per_cpu"),
    threads: config.get("default_resources").get("threads"),
    log:
        "results/mosdepth/{sample}_{type}_coverage_per_amplicon.log",
    benchmark:
        repeat("results/mosdepth/{sample}_{type}_coverage_per_amplicon.benchmark.tsv",
            config.get("mosdepth_merge", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("mosdepth_merge", {}).get("container", config["default_container"])
    message:
        "{rule}: Create merged report for mosdepth"
    script:
        "../scripts/mosdepth_merge.py"


rule mosdepth_overlap_timestep:
    input:
        bamdir = os.path.join(config["runfolder"], "{sample}", config["runid"], "bam_pass"),
        amplibed = [
            f"{config.get('bed_files')}/{target}.bed"
            for target in config.get("amplicons") + config.get("extra_regions")
        ],
    output:
        outdir = temp(directory("results/mosdepth/timestep/{sample}")),
    resources:
        partition=config.get("mosdepth", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("mosdepth", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("mosdepth", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("mosdepth", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mosdepth", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("mosdepth", {}).get("threads", config["default_resources"]["threads"]),
    log:
        "results/mosdepth/timestep/{sample}.log",
    benchmark:
        repeat("results/mosdepth/timestep/{sample}.tsv",
            config.get("mosdepth", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("mosdepth", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Compute coverage with mosdepth for each amplicon.
        """
    script:
        "../scripts/process_timestep_data.py"

rule mosdepth_merge_timestep:
    input:
        indir="results/mosdepth/timestep/{sample}",
    output:
        outdir=temp(directory("results/mosdepth/timestep_coverage/{sample}")),
    resources:
        partition=config.get("default_resources").get("partition"),
        time=config.get("default_resources").get("time"),
        threads=config.get("default_resources").get("threads"),
        mem_mb=config.get("default_resources").get("mem_mb"),
        mem_per_cpu=config.get("default_resources").get("mem_per_cpu"),
    threads: config.get("default_resources").get("threads")
    container:
        config.get("mosdepth_merge", {}).get("container", config["default_container"])
    log:
        "results/mosdepth/timestep/{sample}/timestep_coverage.log",
    benchmark:
        repeat("results/mosdepth/timestep/{sample}/timestep_coverage.benchmark.tsv",
            config.get("mosdepth_merge_timestep", {}).get("benchmark_repeats", 1),
        ),
    message:
        "{rule}: Create merged report for mosdepth"
    script:
        "../scripts/mosdepth_merge_timestep.py"


# rule copy_mosdepth_merge_timestep:
#     input:
#         unpack(get_bam_pass_sample),
#     output:
#         outdir=temp(directory("results/mosdepth/timestep_coverage")),
#         # csv=temp("results/mosdepth/timestep_coverage/timestep{nbatch}_coverage_per_amplicon.csv"),
#     wildcard_constraints:
#         fname=r"[A-Z]{3}\d{3,5}_pass_([a-z0-9]+_)+",
#         nbatch=r"\d+",
#     resources:
#         partition=config["default_resources"]["partition"],
#         time=config["default_resources"]["time"],
#         threads=config["default_resources"]["threads"],
#         mem_mb=config["default_resources"]["mem_mb"],
#         mem_per_cpu=config["default_resources"]["mem_per_cpu"]
#     threads: config["default_resources"]["threads"]
#     log:
#         "results/mosdepth/timestep_coverage/copy_timestep_coverage.log",
#     benchmark:
#         repeat("results/mosdepth/timestep_coverage/copy_timestep_coverage.benchmark.tsv",
#             config.get("copy_mosdepth_merge_timestep", {}).get("benchmark_repeats", 1),
#         )
#     container:
#         config["default_container"]
#     message:
#         """
#         {rule}: Save all timestepped coverage with mosdepth for each amplicon.
#         """
#     shell:
#         """
#         mkdir -p {output.outdir} 2> {log}
#         for fcsv in {input}
#         do
#             cp $fcsv {output.outdir}/$( basename $fcsv ) 2>> {log}
#         done
#         """


rule plot_yield_timestep:
    input:
        indir="results/mosdepth/timestep_coverage/{sample}",
    output:
        csv=temp("results/mosdepth/timestep_coverage_images/{sample}_{type}_cumsum_coverage_per_amplicon.csv"),
        png=temp("results/mosdepth/timestep_coverage_images/{sample}_{type}_cumsum_coverage_per_amplicon.png"),
    resources:
        partition=config["default_resources"]["partition"],
        time=config["default_resources"]["time"],
        threads=config["default_resources"]["threads"],
        mem_mb=config["default_resources"]["mem_mb"],
        mem_per_cpu=config["default_resources"]["mem_per_cpu"],
    threads: config["default_resources"]["threads"]
    container:
        config.get("plot_yield_timestep", {}).get("container", config["default_container"])
    log:
        "results/mosdepth/timestep_coverage_images/{sample}_{type}_cumsum_coverage_per_amplicon.log",
    benchmark:
        repeat("results/mosdepth/timestep_coverage_images/{sample}_{type}_cumsum_coverage_per_amplicon.benchmark.tsv",
            config.get("plot_yield_timestep", {}).get("benchmark_repeats", 1),
        )
    message:
        """
        {rule}: Plot sequencing output per amplicon over time.
        """
    script:
        "../scripts/seq_yield_timestep.py"

rule sequali:
    input:
        fastgz1="prealignment/filtlong/{sample}_{type}_reads.ont_adapt_trim.filtered.fastq.gz",
        fastgz2= "prealignment/filtlong/{sample}_{type}_reads.ont_adapt_trim.filtered.out.fastq.gz",
    output:
        html1=temp("results/sequali/{sample}_{type}_reads.ont_adapt_trim.filtered.fastq.gz.html"),
        json1=temp("results/sequali/{sample}_{type}_reads.ont_adapt_trim.filtered.fastq.gz.json"),
        html2= temp("results/sequali/{sample}_{type}_reads.ont_adapt_trim.filtered.out.fastq.gz.html"),
        json2= temp("results/sequali/{sample}_{type}_reads.ont_adapt_trim.filtered.out.fastq.gz.json"),
    resources:
        partition=config.get("sequali", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("sequali", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("sequali", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("sequali", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sequali", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("sequali", {}).get("threads", config["default_resources"]["threads"])
    log:
        "results/sequali/{sample}_{type}_sequali.log",
    benchmark:
        repeat("results/sequali/{sample}_{type}_sequali.benchmark.tsv",
            config.get("sequali", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("sequali", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Parse the (out)filtered reads and create quality report with sequali.
        """
    shell:
        """
        sequali --html {output.html1} --json {output.json1} {input.fastgz1} 2> {log}
        sequali --html {output.html2} --json {output.json2} {input.fastgz2} 2>> {log}
        """

rule yield_per_pool:
    input:
        csv="results/mosdepth/{sample}_{type}_coverage_per_amplicon.csv",
    output:
        csv=temp("results/mosdepth/{sample}_{type}_yield_pool_{pooln}.csv"),
    resources:
        partition=config.get("default_resources").get("partition"),
        time=config.get("default_resources").get("time"),
        threads=config.get("default_resources").get("threads"),
        mem_mb=config.get("default_resources").get("mem_mb"),
        mem_per_cpu=config.get("default_resources").get("mem_per_cpu"),
    threads: config.get("default_resources").get("threads")
    container:
        config.get("yield_per_pool", {}).get("container", config["default_container"])
    log:
        "results/mosdepth/{sample}_{type}_yield_pool_{pooln}.log",
    benchmark:
        repeat("results/mosdepth/{sample}_{type}_yield_pool_{pooln}.benchmark.tsv",
            config.get("yield_per_pool", {}).get("benchmark_repeats", 1),
        )
    message:
        "{rule}: Calculate number of reads per pool"
    script:
        "../scripts/yield_per_pool.py"

# Picard HsMetrics requires interval file
rule bed_to_interval_list:
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
        repeat("results/qc/picard/BedToIntervalList.benchmark.tsv",
            config.get("bed_to_interval_list", {}).get("benchmark_repeats", 1),
        )
    container:
        config.get("picard_bed_to_interval_list", {}).get("container", config["default_container"])
    wrapper:
        "v5.0.1/bio/picard/bedtointervallist"
