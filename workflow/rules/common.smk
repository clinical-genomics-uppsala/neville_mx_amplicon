__author__ = "Camille Clouard"
__copyright__ = "Copyright 2024, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"

import itertools
import numpy as np
import pathlib
import pandas as pd
import yaml
import sys
import os
import re
from datetime import datetime
from snakemake.utils import validate
from snakemake.utils import min_version
from snakemake.logging import logger

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *

from hydra_genetics.utils.misc import export_config_as_file, get_module_snakefile
from hydra_genetics.utils.software_versions import add_version_files_to_multiqc
from hydra_genetics.utils.software_versions import add_software_version_to_config
from hydra_genetics.utils.software_versions import export_pipeline_version_as_file
from hydra_genetics.utils.software_versions import export_software_version_as_file
from hydra_genetics.utils.software_versions import get_pipeline_version
from hydra_genetics.utils.software_versions import touch_pipeline_version_file_name
from hydra_genetics.utils.software_versions import touch_software_version_file
from hydra_genetics.utils.software_versions import use_container

min_version("7.8.0")

logger.info(f"\n{workflow.snakefile} is being parsed")

### Set and validate config file

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")

try:
    validate(config, schema="../schemas/config.schema.yaml")
except WorkflowError as we:
    # Probably a validation error, but the original exception in lost in
    # snakemake. Pull out the most relevant information instead of a potentially
    # *very* long error message.
    if not we.args[0].lower().startswith("error validating config file"):
        raise
    error_msg = "\n".join(we.args[0].splitlines()[:2])
    parent_rule_ = we.args[0].splitlines()[3].split()[-1]
    if parent_rule_ == "schema:":
        sys.exit(error_msg)
    else:
        schema_hiearachy = parent_rule_.split()[-1]
        schema_section = ".".join(re.findall(r"\['([^']+)'\]", schema_hiearachy)[1::2])
        sys.exit(f"{error_msg} in {schema_section}")

pipeline_version = get_pipeline_version(workflow, pipeline_name="neville_mx_amplicon")



### Read and validate resources file

config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

try:
    samples = pd.read_table(config["samples_run"], dtype=str)["sample"].set_index("sample", drop=False)
except AttributeError:
    samples = pd.read_table(config["samples_run"], dtype=str)["sample"].to_frame().set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")


### Read and validate units file

units = (
    pandas.read_table(config["units_run"], dtype=str)
    .set_index(
        [
            "sample",
            "type",
            "platform",
            "machine",
            "processing_unit",
            "run_id",
            "barcode",
            "methylation",
            "basecalling_model",
            "bam",
        ],
        drop=False,
    )
    .sort_index()
)

validate(units, schema="../schemas/units.schema.yaml")

### Read and validate output file

with open(config["output"]) as output:
    if config["output"].endswith("json"):
        output_spec = json.load(output)
    elif config["output"].endswith("yaml") or config["output"].endswith("yml"):
        output_spec = yaml.safe_load(output.read())

validate(output_spec, schema="../schemas/output_files.schema.yaml")


### Set wildcard constraints
wildcard_constraints:
    experiment=config["batchid"],
    sample="|".join(samples.index),
    type="N|T",
    report="amplicons",
    target="|".join(config.get("amplicons") + config.get("extra_regions")).replace('+', '\+'),  # escape the '+' which has a specific meaning in regex


logger.info("\nTargets: " + "|".join(config.get("amplicons") + config.get("extra_regions")))

### Define functions to be used in the workflow


def read_bam_pass_names(*args):
    bamdir = os.path.join(*args)
    names = []
    batches = []
    for bfile in os.listdir(bamdir):
        if bfile.endswith(".bam"):
            name = '_'.join(bfile.split('_')[:-1])
            if name not in names:
                names.append(name)
            batches.append(bfile.split('_')[-1].replace(".bam", ""))
    return names, batches


def get_bam_pass_sample():
    expr = lambda wildcards: expand(
        "results/mosdepth/timestep/{{fname}}_{{nbatch}}/{target}.mosdepth.summary.txt",
        fname=read_bam_pass_names(config["runfolder"], f"{wildcards.sample}", config["runid"], "bam_pass")[0],
        nbatch=read_bam_pass_names(config["runfolder"], f"{wildcards.sample}", config["runid"], "bam_pass")[1],
        target=config.get("amplicons") + config.get("extra_regions"),
    )
    return expr


def compile_output_file_list(wildcards):
    outdir = pathlib.Path(output_spec.get("directory", "./"))
    output_files = []

    for f in output_spec["files"]:
        # Please remember to add any additional values down below
        # that the output strings should be formatted with.
        outputpaths = set(
            [
                f["output"].format(
                    experiment=experiment, sample=sample, type=unit_type, target=target, report=report, caller=caller
                )
                for experiment in [config["batchid"]]
                for sample in get_samples(samples)
                for unit_type in get_unit_types(units, sample)
                for target in config.get("amplicons") + config.get("extra_regions")
                for report in config["multiqc"]["reports"]
                for caller in config.get("callers", ["clairs_to", "deepsomatic", "vardict"])
            ]
        )

        for op in outputpaths:
            output_files.append(outdir / Path(op))

    return output_files


def generate_copy_rules(output_spec):
    output_directory = pathlib.Path(output_spec.get("directory", "./"))
    rulestrings = []

    for f in output_spec["files"]:
        if f["input"] is None:
            continue

        rule_name = "_copy_{}".format("_".join(re.split(r"\W{1,}", f["name"].strip().lower())))
        input_file = pathlib.Path(f["input"])
        if config.get("subset_vcf_bed"):
            input_str = str(input_file)
            if input_str.endswith(".vcf"):
                input_file = pathlib.Path(input_str[:-4] + ".subset.vcf")
            elif input_str.endswith(".vcf.gz"):
                input_file = pathlib.Path(input_str[:-7] + ".subset.vcf.gz")
            elif input_str.endswith(".vcf.gz.tbi"):
                input_file = pathlib.Path(input_str[:-11] + ".subset.vcf.gz.tbi")
        output_file = output_directory / pathlib.Path(f["output"])
        if pathlib.Path(f["output"]).suffix == '':  # no file extension, hence a directory
            _ = f'@workflow.output(directory("{output_file}"))\n'
        else:
            _ = f'@workflow.output("{output_file}")\n'

        mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
        mem_per_cpu = config.get("_copy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])
        partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
        threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
        time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
        copy_container = config.get("_copy", {}).get("container", config["default_container"])

        rule_code = "\n".join(
            [
                f'@workflow.rule(name="{rule_name}")',
                f'@workflow.input("{input_file}")',
                f'@workflow.output(directory("{output_file}"))'
                if pathlib.Path(f["output"]).suffix == ''
                else f'@workflow.output("{output_file}")',
                f'@workflow.log("logs/{rule_name}_{output_file.name}.log")',
                f'@workflow.container("{copy_container}")',
                f'@workflow.resources(time="{time}", threads={threads}, mem_mb="{mem_mb}", '
                f'mem_per_cpu={mem_per_cpu}, partition="{partition}")',
                f'@workflow.shellcmd("{copy_container}")',  # replace with '@workflow.shellcmd("cp -r {input} {output}")\n\n' ?
                "@workflow.run\n",
                f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, "
                "log, version, rule, conda_env, container_img, singularity_args, use_singularity, "
                "env_modules, bench_record, jobid, is_shell, bench_iteration, cleanup_scripts, "
                "shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
                "__is_snakemake_rule_func=True):",
                '\tshell("(cp -r --preserve=timestamps {input[0]} {output[0]}) &> {log}", bench_record=bench_record, '
                "bench_iteration=bench_iteration)\n\n",
            ]
        )

        rulestrings.append(rule_code)

    exec(compile("\n".join(rulestrings), "copy_result_files", "exec"), workflow.globals)


generate_copy_rules(output_spec)


rule software_versions:
    output:
        software_versions="results/versions/software/softwares_mqc_versions.yaml",
        pipeline_versions="results/versions/software/pipeline_versions_mqc_versions.yaml"
    log:
        "results/versions/software/softwares_mqc_versions.log"
    run:
        import yaml
        import logging
        from datetime import datetime

        logger = logging.getLogger(__name__)

        # 1. Extract pipeline version
        pipeline_info = get_pipeline_version(workflow, pipeline_name="neville_mx_amplicon")
        
        # 2. Extract software versions from containers using hydra-genetics utility
        software_info = {}
        if use_container(workflow):
            try:
                _, software_info = add_software_version_to_config(config, workflow, False)
            except Exception as e:
                logger.warning(f"Could not extract software versions with add_software_version_to_config: {e}")

        # 3. Fallback extraction logic for all containers defined in config
        def get_all_containers(d, containers=None):
            if containers is None:
                containers = {}
            for k, v in d.items():
                if isinstance(v, dict):
                    get_all_containers(v, containers)
                elif k in ["container", "default_container"]:
                    val = str(v)
                    if "/" in val:
                        parts = val.split("/")[-1].split(":")
                        name = parts[0]
                        version = parts[1] if len(parts) > 1 else "latest"
                    else:
                        name = val
                        version = "unknown"
                    containers[name] = version
            return containers

        fallback_versions = get_all_containers(config)

        # 4. Consolidate software versions into a clean flat dictionary
        final_software_versions = {}
        if software_info:
            for container_key, versions in software_info.items():
                for tool, ver in versions.items():
                    if tool != "NOTE":
                        final_software_versions[tool] = ver

        # Fill in missing versions from container URIs
        for tool, ver in fallback_versions.items():
            if tool not in final_software_versions:
                final_software_versions[tool] = ver

        # 5. Write the consolidated software versions file
        with open(output.software_versions, "w") as f:
            yaml.dump(final_software_versions, f)

        # 6. Write the pipeline version file
        flat_pipeline_versions = {}
        for p_name, p_data in pipeline_info.items():
            flat_pipeline_versions[p_name] = p_data.get("version", "unknown")
            
        with open(output.pipeline_versions, "w") as f:
            yaml.dump(flat_pipeline_versions, f)

        # Write log
        with open(log[0], "w") as f:
            f.write(f"Software versions extracted successfully at {datetime.now()}\n")
            f.write(f"Pipeline info: {pipeline_info}\n")
            f.write(f"Software info: {software_info}\n")
            f.write(f"Fallback versions: {fallback_versions}\n")
            f.write(f"Final consolidated versions: {final_software_versions}\n")
