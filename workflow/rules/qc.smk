__author__ = "Camille Clouard"
__copyright__ = "Copyright 2024, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"

import os
def read_bam_pass_names(bamdir):
    names = []
    batches = []
    for bfile in os.listdir(bamdir):
        if bfile.endswith(".bam"):
            names.append('_'.join(bfile.split('_')[:-1]))
            batches.append(bfile.split('_')[-1].replace(".bam", ""))
    return names, batches

bam_pass, batches = read_bam_pass_names(os.path.join(config["runfolder"], "bam_pass"))

wildcard_constraints:
    fname="[A-Z]{3}\d{3}_pass_[a-z0-9]+_[a-z0-9]+",
    nbatch="\d+"

rule pycoqc:
    input:
        seq_run_dir = config.get("runfolder")
    output:
        html = temp("results/pycoqc/{sample}_{type}_report_sequencing_summary.html"),
        json = temp("results/pycoqc/{sample}_{type}_report_sequencing_summary.json"),
        txt = temp("results/pycoqc/{sample}_{type}_report_sequencing_summary.txt")
    resources:
        partition=config.get("pycoqc",{}).get("partition",config["default_resources"]["partition"]),
        time=config.get("pycoqc",{}).get("time",config["default_resources"]["time"]),
        threads=config.get("pycoqc",{}).get("threads",config["default_resources"]["threads"]),
        mem_mb=config.get("pycoqc",{}).get("mem_mb",config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pycoqc",{}).get("mem_per_cpu",config["default_resources"]["mem_per_cpu"]),
    container:
        config.get("pycoqc",{}).get("container",config["default_container"])
    message:
        """
        {rule}: Report graphically for the sequencing run.
        """
    shell:
        """
        summary=$( ls {input.seq_run_dir}/sequencing_summary*.txt ) 
        cp $summary {output.txt}
        pycoQC -f {output.txt} --html_outfile {output.html} --json_outfile {output.json}
        """

rule mosdepth_overlap:
    input:
        bam = "alignment/dorado_align/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam",
        bamidx = "alignment/dorado_align/{sample}_{type}_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam.bai",
        amplibed = os.path.join(config.get("bed_files"), "{target}.bed")
    output:
        bed=temp("results/mosdepth/{sample}_{type}_{target}.regions.bed.gz"),
        csi=temp("results/mosdepth/{sample}_{type}_{target}.regions.bed.gz.csi"),
        glob=temp("results/mosdepth/{sample}_{type}_{target}.mosdepth.global.dist.txt"),
        region=temp("results/mosdepth/{sample}_{type}_{target}.mosdepth.region.dist.txt"),
        summary=temp("results/mosdepth/{sample}_{type}_{target}.mosdepth.summary.txt")
    params:
        prefix_out = 'results/mosdepth',
        threads = 20
    resources:
        partition = config.get("mosdepth", {}).get("partition", config["default_resources"]["partition"]),
        time = config.get("mosdepth", {}).get("time", config["default_resources"]["time"]),
        threads = config.get("mosdepth", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb = config.get("mosdepth", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu = config.get("mosdepth", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    container:
        config.get("mosdepth", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Compute coverage with mosdepth for each amplicon.
        """
    shell:
        """
        chrom=$( cat {input.amplibed} | cut -d$'\t' -f1 )
        mosdepth -t {params.threads} -c $chrom -b {input.amplibed}  {params.prefix_out}/{wildcards.sample}_{wildcards.type}_{wildcards.target} {input.bam}
        """

rule mosdepth_merge:
    input:
        expand("results/mosdepth/{{sample}}_{{type}}_{target}.mosdepth.summary.txt",
            target=config.get("amplicons") + config.get("extra_regions"))
    output:
        csv=temp("results/mosdepth/{sample}_{type}_coverage_per_amplicon.csv")
    resources:
        partition = config.get("default_resources").get("partition"),
        time = config.get("default_resources").get("time"),
        threads = config.get("default_resources").get("threads"),
        mem_mb = config.get("default_resources").get("mem_mb"),
        mem_per_cpu = config.get("default_resources").get("mem_per_cpu"),
    message:
        "{rule}: Create merged report for mosdepth"
    run:
        import os
        import pandas as pd

        cols = ["length", "bases", "mean", "target"]
        df = pd.concat([pd.read_csv(summary, sep='\t')
                        .assign(target=os.path.basename(summary).replace(".mosdepth.summary.txt", "")) \
                        for summary in list(input)
                        ]).loc[3,cols] # lines for "total_region" have index=3 in the dataframe
        df.set_index("target", inplace=True)
        df.to_csv(output.csv, index=True)

rule mosdepth_overlap_timestep:
    input:
        bam = os.path.join(config["runfolder"],
            "bam_pass/{fname}_{nbatch}.bam"),
        bamidx = os.path.join(config["runfolder"],
            "bam_pass/{fname}_{nbatch}.bam.bai"),
        amplibed = os.path.join(config.get("bed_files"), "{target}.bed")
    output:
        bed=temp("results/mosdepth/timestep/{fname}_{nbatch}/{target}.regions.bed.gz"),
        csi=temp("results/mosdepth/timestep/{fname}_{nbatch}/{target}.regions.bed.gz.csi"),
        glob=temp("results/mosdepth/timestep/{fname}_{nbatch}/{target}.mosdepth.global.dist.txt"),
        region=temp("results/mosdepth/timestep/{fname}_{nbatch}/{target}.mosdepth.region.dist.txt"),
        summary=temp("results/mosdepth/timestep/{fname}_{nbatch}/{target}.mosdepth.summary.txt")
    params:
        prefix_out = 'results/mosdepth/timestep',
    resources:
        partition = config.get("mosdepth", {}).get("partition", config["default_resources"]["partition"]),
        time = config.get("mosdepth", {}).get("time", config["default_resources"]["time"]),
        threads = config.get("mosdepth", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb = config.get("mosdepth", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu = config.get("mosdepth", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    container:
        config.get("mosdepth", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Compute coverage with mosdepth for each amplicon.
        """
    shell:
        """
        chrom=$( cat {input.amplibed} | cut -d$'\t' -f1 )
        mosdepth -t {resources.threads} -c $chrom -b {input.amplibed}  {params.prefix_out}/{wildcards.fname}_{wildcards.nbatch}/{wildcards.target} {input.bam}
        """

rule mosdepth_merge_timestep:
    input:
        expand("results/mosdepth/timestep/{{fname}}_{{nbatch}}/{target}.mosdepth.summary.txt",
            fname=bam_pass,
            nbatch=batches,
            target=config.get("amplicons") + config.get("extra_regions"))
    output:
        csv=temp("results/mosdepth/timestep/{fname}_{nbatch}/timestep{nbatch}_coverage_per_amplicon.csv"),
    resources:
        partition = config.get("default_resources").get("partition"),
        time = config.get("default_resources").get("time"),
        threads = config.get("default_resources").get("threads"),
        mem_mb = config.get("default_resources").get("mem_mb"),
        mem_per_cpu = config.get("default_resources").get("mem_per_cpu"),
    message:
        "{rule}: Create merged report for mosdepth"
    run:
        import os
        import pandas as pd

        cols = ["length", "bases", "mean", "target"]
        df = pd.concat([pd.read_csv(summary, sep='\t')
                        .assign(target=os.path.basename(summary).replace(".mosdepth.summary.txt", "")) \
                        for summary in list(input)
                        ])
        print(df)
        try:
            df = df[df["chrom"] == "total_region"].loc[:, cols]  # .loc[3,cols] # lines for "total_region" have index=3 in the dataframe
        except: # this might not be necessary
            df = pd.DataFrame([0.0, 0.0, 0.0, 0.0], columns=cols)
        print(df)
        try:
            df.set_index("target", inplace=True)
        except AttributeError:
            df.rename("target", inplace=True)
        df.to_csv(output.csv, index=True)


rule copy_mosdepth_merge_timestep:
    input:
        expand("results/mosdepth/timestep/{fname}_{nbatch}/timestep{nbatch}_coverage_per_amplicon.csv",
            fname=bam_pass,
            nbatch=batches,
        )
    output:
        outdir=temp(directory("results/mosdepth/timestep_coverage")),
        # csv=temp("results/mosdepth/timestep_coverage/timestep{nbatch}_coverage_per_amplicon.csv")
    resources:
        partition=config["default_resources"]["partition"],
        time=config["default_resources"]["time"],
        threads=config["default_resources"]["threads"],
        mem_mb=config["default_resources"]["mem_mb"],
        mem_per_cpu=config["default_resources"]["mem_per_cpu"]
    # container:
    #     config["default_container"]
    message:
        """
        {rule}: Save all timestepped coverage with mosdepth for each amplicon.
        """
    shell:
        """
        mkdir -p {output.outdir}
        for fcsv in {input}
        do
            cp $fcsv {output.outdir}/$( basename $fcsv )
        done
        """


rule plot_yield_timestep:
    input:
        indir="results/mosdepth/timestep_coverage"
    output:
        csv=temp("results/mosdepth/timestep_coverage_images/{sample}_{type}_cumsum_coverage_per_amplicon.csv"),
        png=temp("results/mosdepth/timestep_coverage_images/{sample}_{type}_cumsum_coverage_per_amplicon.png")
    resources:
        partition=config["default_resources"]["partition"],
        time=config["default_resources"]["time"],
        threads=config["default_resources"]["threads"],
        mem_mb=config["default_resources"]["mem_mb"],
        mem_per_cpu=config["default_resources"]["mem_per_cpu"]
    # container:
    #     config["default_container"]
    message:
        """
        {rule}: Plot sequencing output per amplicon over time.
        """
    script:
        "../scripts/seq_yield_timestep.py"

rule sequali:
    input:
        fastgz1 = "prealignment/filtlong/{sample}_{type}_reads.ont_adapt_trim.filtered.fastq.gz",
        fastgz2= "prealignment/filtlong/{sample}_{type}_reads.ont_adapt_trim.filtered.out.fastq.gz"
    output:
        html1 = temp("results/sequali/{sample}_{type}_reads.ont_adapt_trim.filtered.fastq.gz.html"),
        json1 = temp("results/sequali/{sample}_{type}_reads.ont_adapt_trim.filtered.fastq.gz.json"),
        html2= temp("results/sequali/{sample}_{type}_reads.ont_adapt_trim.filtered.out.fastq.gz.html"),
        json2= temp("results/sequali/{sample}_{type}_reads.ont_adapt_trim.filtered.out.fastq.gz.json"),
    resources:
        partition=config.get("sequali",{}).get("partition",config["default_resources"]["partition"]),
        time=config.get("sequali",{}).get("time",config["default_resources"]["time"]),
        threads=config.get("sequali",{}).get("threads",config["default_resources"]["threads"]),
        mem_mb=config.get("sequali",{}).get("mem_mb",config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sequali",{}).get("mem_per_cpu",config["default_resources"]["mem_per_cpu"]),
    # container:
    #     config.get("sequali", {}).get("container", config["default_container"])
    message:
        """
        {rule}: Parse the (out)filtered reads and create quality report with sequali.
        """
    shell:
        """
        sequali --html {output.html1} --json {output.json1} {input.fastgz1}
        sequali --html {output.html2} --json {output.json2} {input.fastgz2}
        """

rule yield_per_pool:
    input:
        csv = "results/mosdepth/{sample}_{type}_coverage_per_amplicon.csv"
    output:
        csv = temp("results/mosdepth/{sample}_{type}_yield_pool_{pooln}.csv")
    resources:
        partition=config.get("default_resources").get("partition"),
        time=config.get("default_resources").get("time"),
        threads=config.get("default_resources").get("threads"),
        mem_mb=config.get("default_resources").get("mem_mb"),
        mem_per_cpu=config.get("default_resources").get("mem_per_cpu"),
    message:
        "{rule}: Calculate number of reads per pool"
    run:
        import pandas as pd

        cols =["pool","target","reads_counts"]
        df = pd.read_csv(input.csv, sep=',', index_col="target")
        counts = dict([(target, round(df.loc[target, "mean"])) \
                       for target in df.index \
                       if target.replace(f"{wildcards.sample}_{wildcards.type}_", "") in config.get("pools")[int(wildcards.pooln)]  # dirty fix to not include the sample's name in the target's label
                       ])
        d_and_j = [target for target in counts.keys() if target.find("+J") >= 0]
        d_only = [(target, round(df.loc[target, "mean"])) \
                  for target in df.index \
                  if (target.find("TP53_D") >= 0 and target.find("_only") >= 0)
                  ]
        if d_and_j:
            counts["TP53_JX"] = counts[d_and_j[0]] - d_only[0][1]
            del counts[d_and_j[0]]
        counts["total"] = sum(counts.values())
        df_p = pd.DataFrame.from_dict(counts,orient="index",columns=["reads_counts"])
        df_p.index.name = f"Pool {wildcards.pooln}"
        df_p = df_p.assign(pct_reads = round(df_p["reads_counts"] * 100 / counts["total"], ndigits=1))
        df_p.to_csv(output.csv,index=True)

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
        # optional parameters
        extra="--SORT true",  # sort output interval list before writing
    resources:
        partition=config.get("default_resources").get("partition"),
        time=config.get("default_resources").get("time"),
        threads=config.get("default_resources").get("threads"),
        mem_mb=config.get("default_resources").get("mem_mb"),
        mem_per_cpu=config.get("default_resources").get("mem_per_cpu"),
    container:
        config.get("picard_bed_to_interval_list", {}).get("container", config["default_container"])
    wrapper:
        "v5.0.1/bio/picard/bedtointervallist"

# rule j3_only_bam_to_bed:
#     input:
#         bam = "alignment/dorado_align/{sample}_{type}_TP53_J3_only_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.bam"
#     output:
#         bed = os.path.join(config.get("bed_files",""), "TP53_J3_only.bed")

# rule samtools_stats:
#     input:
#         bam="{run_date}_{run_id}/data/reads.ont_adapt_trim.filtered.aligned.sorted.bam",
#         bed=os.path.join(config.get("bed_files"), "amplicons.bed")
#     output:
#         "results/{run_date}_{run_id}/samtools/samtools_stats.txt"
#     resources:
#         partition=config.get("samtools",{}).get("partition",config["default_resources"]["partition"]),
#         time=config.get("samtools",{}).get("time",config["default_resources"]["time"]),
#         threads=config.get("samtools",{}).get("threads",config["default_resources"]["threads"]),
#         mem_mb=config.get("samtools",{}).get("mem_mb",config["default_resources"]["mem_mb"]),
#         mem_per_cpu=config.get("samtools",{}).get("mem_per_cpu",config["default_resources"]["mem_per_cpu"]),
#     container:
#         config.get("samtools",{}).get("container",config["default_container"])
#     shell:
#         """
#         samtools stats -t {input.bed} {input.bam} > {output}
#         """

