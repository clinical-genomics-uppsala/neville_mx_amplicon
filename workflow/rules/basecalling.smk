__author__ = "Camille Clouard"
__copyright__ = "Copyright 2024, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"

import os

# Rules to be run on single-sequenced samples
if not config.get("multisample", False):


    rule basecalling_dorado:
        input:
            pod5 = os.path.join(
                config.get("runfolder"),
                config.get("batchid"),
                config.get("runid"),
                config.get("raw_data")
            ),
            ref_data = config.get("ref_data")
        output:
            bam = temp("basecalling/dorado/{sample}_{type}_reads.ont_adapt_trim.bam")
        params:
            dir_models = config.get("dir_models"),  # os.path.abspath(os.pardir),
            dorado_model = config.get("dorado_model"),
            dorado_options = "--device cuda:all --trim adapters"  # --trim all
        resources:
            partition=config.get("dorado_basecalling", {}).get("partition", config["default_resources"]["partition"]),
            time=config.get("dorado_basecalling", {}).get("time", config["default_resources"]["time"]),
            gres=config.get("dorado_basecalling", {}).get("gres"),
            threads=config.get("dorado_basecalling", {}).get("threads", config["default_resources"]["threads"]),
            mem_mb=config.get("dorado_basecalling", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("dorado_basecalling", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
            slurm_extra=config.get("dorado_basecalling", {}).get("slurm_extra"),
        threads: config.get("dorado_basecalling", {}).get("threads", config["default_resources"]["threads"]),
        container:
            config.get("dorado", {}).get("container", config["default_container"])
        log:
            "basecalling/dorado/{sample}_{type}_reads.ont_adapt_trim.bam.log"
        benchmark:
            repeat(
                "basecalling/dorado/{sample}_{type}_reads.ont_adapt_trim.bam.benchmark.tsv",
                config.get("dorado_basecalling", {}).get("benchmark_repeats", 1)
            )
        message:
            "{rule}: Basecalling with dorado from POD5 files. ONT adapters will be trimmed."
        shell:
            """
            echo "Dorado executed from $( which dorado )"
    
            echo "Executing dorado basecalling in {input.pod5} with options '{params.dorado_options}'"
            echo "and model {params.dorado_model}"
    
            dorado basecaller {params.dir_models}/{params.dorado_model} {params.dorado_options} {input.pod5}/ > {output.bam} 2>> {log}
            """


    rule duplex_basecalling_dorado:
        input:
            pod5=os.path.join(
                config.get("runfolder"),
                config.get("batchid"),
                config.get("runid"),
                config.get("raw_data")
            ),
            ref_data=config.get("ref_data"),
        output:
            bam=temp("basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.bam")
        params:
            dir_models=config.get("dir_models"),
            dorado_model=config.get("dorado_model"),
            dorado_options="--device cuda:all",
            trim_options="--no-trim-primers --sequencing-kit SQK-LSK114"
        resources:
            partition=config.get("dorado_basecalling", {}).get("partition", config["default_resources"]["partition"]),
            time=config.get("dorado_basecalling", {}).get("time", config["default_resources"]["time"]),
            gres=config.get("dorado_basecalling", {}).get("gres"),
            threads=config.get("dorado_basecalling", {}).get("threads", config["default_resources"]["threads"]),
            mem_mb=config.get("dorado_basecalling", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("dorado_basecalling", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
            slurm_extra=config.get("dorado_basecalling", {}).get("slurm_extra"),
        threads: config.get("dorado_basecalling", {}).get("threads", config["default_resources"]["threads"]),
        benchmark:
            repeat(
                "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.bam.benchmark.tsv",
                config.get("dorado_basecalling", {}).get("benchmark_repeats", 1)
            )
        container:
            config.get("dorado", {}).get("container", config["default_container"])
        log:
            "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.bam.log"
        message:
            "{rule}: Duplex basecalling with dorado from POD5 files. ONT adapters will be trimmed."
        shell:
            """
            echo "Dorado executed from $( which dorado )" > {log}
    
            echo "Executing dorado duplex basecalling in {input.pod5} with options '{params.dorado_options}'" >> {log}
            echo "and model {params.dorado_model}" >> {log}
            echo "POD5 files found:"
            ls -la {input.pod5}/ >> {log}
    
            dorado duplex {params.dir_models}/{params.dorado_model} {params.dorado_options} {input.pod5}/ 2>> {log} | dorado trim {params.trim_options} > {output.bam} 2>> {log}
            """


# Rules to be run on multiplexed samples
if config.get("multisample", False):


    rule duplex_basecalling_multisamples_dorado:
        input:
            pod5=os.path.join(
                config.get("runfolder"),
                config.get("batchid"),
                config.get("runid"),
                config.get("raw_data")
            ),
            ref_data=config.get("ref_data"),
        output:
            bam = temp("basecalling/dorado_duplex_multisamples/multi_samples_reads.basecalled.bam")
        params:
            dir_models=config.get("dir_models"),
            dorado_model=config.get("dorado_model"),
            dorado_options="--device cuda:all",
        resources:
            partition=config.get("dorado_basecalling", {}).get("partition", config["default_resources"]["partition"]),
            time=config.get("dorado_basecalling", {}).get("time", config["default_resources"]["time"]),
            gres=config.get("dorado_basecalling", {}).get("gres"),
            threads=config.get("dorado_basecalling", {}).get("threads", config["default_resources"]["threads"]),
            mem_mb=config.get("dorado_basecalling", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("dorado_basecalling", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
            slurm_extra=config.get("dorado_basecalling", {}).get("slurm_extra"),
        threads: config.get("dorado_basecalling", {}).get("threads", config["default_resources"]["threads"]),
        benchmark:
            repeat(
                "basecalling/dorado_duplex_multisamples/multi_samples_reads.basecalled.bam.benchmark.tsv",
                config.get("dorado_basecalling", {}).get("benchmark_repeats", 1)
            )
        container:
            config.get("dorado", {}).get("container", config["default_container"])
        log:
            "basecalling/dorado_duplex_multisamples/multi_samples_reads.basecalled.bam.log"
        message:
            "{rule}: Duplex basecalling with dorado inmultiplexed samples from POD5 files. ONT adapters will NOT be trimmed."
        shell:
            """
            echo "Dorado executed from $( which dorado )" > {log}

            echo "Executing dorado duplex basecalling in {input.pod5} with options '{params.dorado_options}'" >> {log}
            echo "and model {params.dorado_model}" >> {log}
            echo "POD5 files found:"
            ls -la {input.pod5}/ >> {log}

            dorado duplex {params.dir_models}/{params.dorado_model} {params.dorado_options} {input.pod5}/ > {output.bam} 2>> {log}
            """


    rule demux_dorado:
        input:
            bam="basecalling/dorado_duplex_multisamples/multi_samples_reads.basecalled.bam",
        output:
            bamdir=temp(directory("basecalling/dorado_demux/")),
        params:
            dorado_options="--kit-name SQK-NBD114-24",
            samplesheet=config.get("samplesheet"),
        resources:
            partition=config.get("trim_dorado",{}).get("partition",config["default_resources"]["partition"]),
            time=config.get("trim_dorado",{}).get("time",config["default_resources"]["time"]),
            threads=config.get("trim_dorado",{}).get("threads",config["default_resources"]["threads"]),
            mem_mb=config.get("trim_dorado",{}).get("mem_mb",config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("trim_dorado",{}).get("mem_per_cpu",config["default_resources"]["mem_per_cpu"]),
        threads: config.get("trim_dorado",{}).get("threads",config["default_resources"]["threads"]),
        benchmark:
            repeat(
                "basecalling/dorado_demux/output.bam.benchmark.tsv",
                config.get("demux_dorado",{}).get("benchmark_repeats",1)
            )
        container:
            config.get("dorado", {}).get("container", config["default_container"])
        log:
            "basecalling/dorado_demux/output.bam.log"
        message:
            "{rule}: Demultiplexing samples with dorado."
        shell:
            """
            echo "Dorado executed from $( which dorado )" > {log}
            echo "Executing dorado demultiplexing in {input.bam} with sample sheet '{params.samplesheet}'" >> {log}
            dorado demux --sample-sheet {params.samplesheet} --output-dir {output.bamdir} {params.dorado_options} {input.bam} &>> {log}
            """


    rule rename_demux_bam:
        input:
            bamdir="basecalling/dorado_demux/",
        output:
            bam_renamed=temp("basecalling/rename_demux_bam/{sample}_{type}_reads.basecalled.bam"),
        resources:
            partition=config.get("rename_demux_bam", {}).get("partition", config["default_resources"]["partition"]),
            time=config.get("rename_demux_bam", {}).get("time", config["default_resources"]["time"]),
            threads=config.get("rename_demux_bam", {}).get("threads", config["default_resources"]["threads"]),
            mem_mb=config.get("rename_demux_bam", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("rename_demux_bam", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        threads: config.get("rename_demux_bam", {}).get("threads", config["default_resources"]["threads"]),
        benchmark:
            repeat(
                "basecalling/rename_demux_bam/{sample}_{type}_reads.basecalled.bam.benchmark.tsv",
                config.get("rename_demux_bam", {}).get("benchmark_repeats", 1)
            )
        container:
            config.get("rename_demux_bam", {}).get("container", config["default_container"])
        log:
            "basecalling/rename_demux_bam/{sample}_{type}_reads.basecalled.bam.log"
        message:
            "{rule}: Renaming demultiplexed BAM files to include sample name and read type."
        shell:
            """
            outdir=$(dirname {output.bam_renamed})
            bams=$(ls {input.bamdir}/*.bam)
            echo $bams > {log}
            for bam in $bams; do
                filename=$(basename -- "$bam")
                sample=$(echo $filename | cut -d'_' -f2 | cut -d'.' -f1)
                if [[ $sample == {wildcards.sample} ]]; then
                    echo "Renaming BAM file for $sample and writing to $outdir" &>> {log} 
                    cp $bam $outdir/${{sample}}_{wildcards.type}_reads.basecalled.bam &>> {log}
                fi
            done
            """


    rule trim_dorado:
        input:
            bam="basecalling/rename_demux_bam/{sample}_{type}_reads.basecalled.bam",
        output:
            bam=temp("basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.bam"),
        params:
            dorado_options="--sequencing-kit SQK-NBD114-24",
        resources:
            partition=config.get("trim_dorado", {}).get("partition", config["default_resources"]["partition"]),
            time=config.get("trim_dorado", {}).get("time", config["default_resources"]["time"]),
            threads=config.get("trim_dorado", {}).get("threads", config["default_resources"]["threads"]),
            mem_mb=config.get("trim_dorado", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("trim_dorado", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        threads: config.get("trim_dorado", {}).get("threads", config["default_resources"]["threads"]),
        benchmark:
            repeat(
                "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.bam.benchmark.tsv",
                config.get("trim_dorado", {}).get("benchmark_repeats", 1)
            )
        container:
            config.get("dorado", {}).get("container", config["default_container"])
        log:
            "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.bam.log"
        message:
            "{rule}: Trimming demultiplexed reads with dorado."
        shell:
            """
            echo "Dorado executed from $( which dorado )" > {log}
    
            echo "Executing dorado trimming in {input.bam} with options '{params.dorado_options}'" >> {log}
    
            dorado trim {input.bam} {params.dorado_options} > {output.bam} 2>> {log}
            """


rule basecalling_bam2fastq:
    input:
        bam="basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.bam",
    output:
        fastq=temp("basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq")
    resources:
        partition=config.get("samtools", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("samtools", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("samtools", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("samtools", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("samtools", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("samtools", {}).get("threads", config["default_resources"]["threads"]),
    benchmark:
        repeat(
            "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq.benchmark.tsv",
            config.get("samtools", {}).get("benchmark_repeats", 1)
        )
    log:
        "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq.log"
    container:
        config.get("samtools", {}).get("container", config["default_container"])
    message:
        "{rule}: Convert unaligned BAM file to FASTQ with samtools"
    shell:
        """
        samtools fastq {input.bam} > {output.fastq} 2> {log}
        """

rule basecalling_compress_fastq:
    input:
        fastq="basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq"
    output:
        fastqgz=temp("basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq.gz")
    resources:
        partition=config.get("samtools", {}).get("partition", config["default_resources"]["partition"]),
        time=config.get("samtools", {}).get("time", config["default_resources"]["time"]),
        threads=config.get("samtools", {}).get("threads", config["default_resources"]["threads"]),
        mem_mb=config.get("samtools", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("samtools", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
    threads: config.get("samtools", {}).get("threads", config["default_resources"]["threads"]),
    benchmark:
        repeat(
            "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq.gz.benchmark.tsv",
            config.get("samtools", {}).get("benchmark_repeats", 1)
        )
    log:
        "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.fastq.gz.log"
    container:
        config.get("basecalling_compress_fastq", {}).get("container", config["default_container"])
    message:
        "{rule}: Compress FASTQ file with basecalled data."
    shell:
        """
        rm -f {output.fastqgz}
        gzip -f {input.fastq} > {output.fastqgz} 2> {log}
        """