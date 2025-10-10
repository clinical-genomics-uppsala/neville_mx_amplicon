__author__ = "Camille Clouard"
__copyright__ = "Copyright 2024, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"

import os

# Rules to be run on single-sequenced samples
if not config.get("multisample", False):


    rule dorado_basecaller:
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
            # dir_models=config.get("dorado_basecaller", {}).get("dirmodels", config.get("dir_models")),
            model=config.get("dorado_basecaller",{}).get("model",""),
            trim=config.get("dorado_basecaller",{}).get("trim",""),
            extra=config.get("dorado_basecaller",{}).get("extra",""),
        resources:
            partition=config.get("dorado_basecaller", {}).get("partition", config["default_resources"]["partition"]),
            time=config.get("dorado_basecaller", {}).get("time", config["default_resources"]["time"]),
            gres=config.get("dorado_basecaller", {}).get("gres"),
            threads=config.get("dorado_basecaller", {}).get("threads", config["default_resources"]["threads"]),
            mem_mb=config.get("dorado_basecaller", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("dorado_basecaller", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
            slurm_extra=config.get("dorado_basecaller", {}).get("slurm_extra"),
        threads: config.get("dorado_basecaller", {}).get("threads", config["default_resources"]["threads"]),
        container:
            config.get("dorado_basecaller", {}).get("container", config["default_container"])
        log:
            "basecalling/dorado/{sample}_{type}_reads.ont_adapt_trim.bam.log"
        benchmark:
            repeat(
                "basecalling/dorado/{sample}_{type}_reads.ont_adapt_trim.bam.benchmark.tsv",
                config.get("dorado_basecaller", {}).get("benchmark_repeats", 1)
            )
        message:
            "{rule}: Basecalling with dorado from POD5 files. ONT adapters will be trimmed."
        shell:
            """
            echo "Dorado executed from $( which dorado )"
    
            echo "Executing dorado basecalling in {input.pod5} with options '{params.trim} {params.extra}'"
            echo "and model {params.model}"
    
            dorado basecaller {params.model} {params.trim} {params.extra} {input.pod5}/ > {output.bam} 2>> {log}
            """


    rule dorado_duplex:
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
            # dir_models=config.get("dorado_duplex", {}).get("dirmodels", config.get("dir_models")),
            model=config.get("dorado_duplex",{}).get("model", ""),
            trim=config.get("dorado_duplex",{}).get("trim",""),
            extra=config.get("dorado_duplex",{}).get("extra",""),
        resources:
            partition=config.get("dorado_basecaller", {}).get("partition", config["default_resources"]["partition"]),
            time=config.get("dorado_basecaller", {}).get("time", config["default_resources"]["time"]),
            gres=config.get("dorado_basecaller", {}).get("gres"),
            threads=config.get("dorado_basecaller", {}).get("threads", config["default_resources"]["threads"]),
            mem_mb=config.get("dorado_basecaller", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("dorado_basecaller", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
            slurm_extra=config.get("dorado_basecaller", {}).get("slurm_extra"),
        threads: config.get("dorado_basecaller", {}).get("threads", config["default_resources"]["threads"]),
        benchmark:
            repeat(
                "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.bam.benchmark.tsv",
                config.get("dorado_basecaller", {}).get("benchmark_repeats", 1)
            )
        container:
            config.get("dorado_duplex", {}).get("container", config["default_container"])
        log:
            "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.bam.log"
        message:
            "{rule}: Duplex basecalling with dorado from POD5 files. ONT adapters will be trimmed."
        shell:
            """
            echo "Dorado executed from $( which dorado )" > {log}
    
            echo "Executing dorado duplex basecalling in {input.pod5} with options '{params.trim} {params.extra}'" >> {log}
            echo "and model {params.model}" >> {log}
            echo "POD5 files found:"
            ls -la {input.pod5}/ >> {log}
    
            dorado duplex {params.model} {params.extra} {input.pod5}/ 2>> {log} | dorado trim {params.trim} > {output.bam} 2>> {log}
            """


# Rules to be run on multiplexed samples
if config.get("multisample", False):


    rule dorado_duplex_multisamples:
        input:
            pod5=os.path.join(
                config.get("runfolder"),
                config.get("batchid"),
                config.get("runid"),
                config.get("raw_data")
            ),
            ref_data=config.get("ref_data"),
        output:
            bam = temp(f"basecalling/dorado_duplex_multisamples/"
                       f"{config['batchid']}/multi_samples_reads.basecalled.bam"),
        params:
            # dir_models=config.get("dorado_duplex_multisamples", {}).get("dirmodels", config.get("dir_models")),
            model=config.get("dorado_duplex_multisamples", {}).get("model", ""),
            trim=config.get("dorado_duplex_multisamples", {}).get("trim", ""),
            extra=config.get("dorado_duplex_multisamples", {}).get("extra", ""),
        resources:
            partition=config.get("dorado_basecaller", {}).get("partition", config["default_resources"]["partition"]),
            time=config.get("dorado_basecaller", {}).get("time", config["default_resources"]["time"]),
            gres=config.get("dorado_basecaller", {}).get("gres"),
            threads=config.get("dorado_basecaller", {}).get("threads", config["default_resources"]["threads"]),
            mem_mb=config.get("dorado_basecaller", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("dorado_basecaller", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
            slurm_extra=config.get("dorado_basecaller", {}).get("slurm_extra"),
        threads: config.get("dorado_basecaller", {}).get("threads", config["default_resources"]["threads"]),
        benchmark:
            repeat(
                "basecalling/dorado_duplex_multisamples/multi_samples_reads.basecalled.bam.benchmark.tsv",
                config.get("dorado_basecaller", {}).get("benchmark_repeats", 1)
            )
        container:
            config.get("dorado_duplex_multisamples", {}).get("container", config["default_container"])
        log:
            "basecalling/dorado_duplex_multisamples/multi_samples_reads.basecalled.bam.log"
        message:
            "{rule}: Duplex basecalling with dorado inmultiplexed samples from POD5 files. ONT adapters will NOT be trimmed."
        shell:
            """
            echo "Dorado executed from $( which dorado )" > {log}
            echo "Downloading model {params.model} if not already present." >> {log}
            dorado download --model {params.model} >> {log}
            echo "Executing dorado duplex basecalling in {input.pod5} with options '{params.trim} {params.extra}'" >> {log}
            echo "and model {params.model}" >> {log}
            echo "POD5 files found:"
            ls -la {input.pod5}/ >> {log}
            dorado duplex {params.model} {params.trim} {params.extra} {input.pod5}/ > {output.bam} 2>> {log}
            rm -rf {params.model}
            """


    rule dorado_demux:
        input:
            bam=f"basecalling/dorado_duplex_multisamples/{config['batchid']}/multi_samples_reads.basecalled.bam",
        output:
            bamdir=temp(directory(os.path.join("basecalling/dorado_demux/", config['batchid']))),
            done=temp(f"basecalling/dorado_demux/{config['batchid']}_demux.done"),
        params:
            extra="--kit-name SQK-NBD114-24",
            samplesheet=config.get("samplesheet"),
        resources:
            partition=config.get("demux_dorado",{}).get("partition",config["default_resources"]["partition"]),
            time=config.get("demux_dorado",{}).get("time",config["default_resources"]["time"]),
            threads=config.get("demux_dorado",{}).get("threads",config["default_resources"]["threads"]),
            mem_mb=config.get("demux_dorado",{}).get("mem_mb",config["default_resources"]["mem_mb"]),
            mem_per_cpu=config.get("demux_dorado",{}).get("mem_per_cpu",config["default_resources"]["mem_per_cpu"]),
        threads: config.get("demux_dorado",{}).get("threads",config["default_resources"]["threads"]),
        benchmark:
            repeat(
                "basecalling/dorado_demux/output.bam.benchmark.tsv",
                config.get("demux_dorado",{}).get("benchmark_repeats",1)
            )
        container:
            config.get("dorado_demux", {}).get("container", config["default_container"])
        log:
            "basecalling/dorado_demux/output.bam.log"
        message:
            "{rule}: Demultiplexing samples with dorado."
        shell:
            """
            echo "Dorado executed from $( which dorado )" > {log}
            echo "Executing dorado demultiplexing in {input.bam} with sample sheet '{params.samplesheet}'" >> {log}
            dorado demux --sample-sheet {params.samplesheet} --output-dir {output.bamdir} {params.extra} {input.bam} &  
            &>> {log}
            process_id=$!
            echo "Waiting for demux to complete... Process PID: $process_id"
            wait $process_id
            touch {output.done}
            """


    rule rename_demux_bam:
        input:
            bamdir=os.path.join("basecalling/dorado_demux/", config['batchid']),
            done=f"basecalling/dorado_demux/{config['batchid']}_demux.done",
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


    rule dorado_trim:
        input:
            bam="basecalling/rename_demux_bam/{sample}_{type}_reads.basecalled.bam",
        output:
            bam=temp("basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.bam"),
        params:
            trim=config.get("dorado_trim", {}).get("trim", ""),
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
            config.get("dorado_trim", {}).get("container", config["default_container"])
        log:
            "basecalling/dorado_duplex/{sample}_{type}_reads.ont_adapt_trim.bam.log"
        message:
            "{rule}: Trimming demultiplexed reads with dorado."
        shell:
            """
            echo "Dorado executed from $( which dorado )" > {log}
    
            echo "Executing dorado trimming in {input.bam} with options '{params.trim}'" >> {log}
    
            dorado trim {input.bam} {params.trim} > {output.bam} 2>> {log}
            """


rule bam2fastq:
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

rule compress_fastq:
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