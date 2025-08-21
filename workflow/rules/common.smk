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

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *

from hydra_genetics.utils.misc import export_config_as_file
from hydra_genetics.utils.software_versions import add_version_files_to_multiqc
from hydra_genetics.utils.software_versions import add_software_version_to_config
from hydra_genetics.utils.software_versions import export_pipeline_version_as_file
from hydra_genetics.utils.software_versions import export_software_version_as_file
from hydra_genetics.utils.software_versions import get_pipeline_version
from hydra_genetics.utils.software_versions import touch_pipeline_version_file_name
from hydra_genetics.utils.software_versions import touch_software_version_file
from hydra_genetics.utils.software_versions import use_container

min_version("7.8.0")

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

date_string = datetime.now().strftime('%Y%m%d--%H-%M-%S')
pipeline_version = get_pipeline_version(workflow, pipeline_name="pipeline_pool_amplicon")
# version_files = touch_pipeline_version_file_name(pipeline_version, date_string=date_string, directory="results/versions/software")
# if use_container(workflow):
#     version_files.append(touch_software_version_file(config, date_string=date_string, directory="results/versions/software"))
# add_version_files_to_multiqc(config, version_files)

onstart:
    export_pipeline_version_as_file(pipeline_version, date_string=date_string, directory="results/versions/software")
    # Make sure that the user have the requested containers to be used
    # if use_container(workflow):
    #     # From the config retrieve all dockers used and parse labels for software versions. Add
    #     # this information to config dict.
    #     update_config, software_info = add_software_version_to_config(config, workflow, False) # fails with Pisces
    #     # Print all softwares used as files. Additional parameters that can be set
    #     # - directory, default value: software_versions
    #     # - file_name_ending, default value: mqc_versions.yaml
    #     # date_string, a string that will be added to the folder name to make it unique (preferably a timestamp)
    #     # export_software_version_as_file(software_info, date_string=date_string, directory="results/versions/software")
    # print config dict as a file. Additional parameters that can be set
    # output_file, default config
    # output_directory, default = None, i.e no folder
    # date_string, a string that will be added to the folder name to make it unique (preferably a timestamp)
    # export_config_as_file(update_config, date_string=date_string, directory="results/versions")

### Read and validate resources file

config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file

try:
    samples = pd.read_table(config["samples_run"], dtype=str)["sample"].set_index("sample", drop=False)
except AttributeError:
    samples = pd.read_table(config["samples_run"],dtype=str)["sample"].to_frame().set_index("sample",drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")
print(samples)

### Read and validate units file

units = (
    pandas.read_table(config["units_run"], dtype=str)
    .set_index(["sample",
                "type",
                "platform",
                "machine",
                "processing_unit",
                "run_id",
                "barcode",
                "methylation",
                "basecalling_model",
                "bam"],
        drop=False)
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
    sample="|".join(samples.index),
    type="N|T",
    report="amplicons",
    target = "|".join(config.get("amplicons") + config.get("extra_regions")).replace('+', '\+')  # escape the '+' which has a specific meaning in regex

print("Targets: ", "|".join(config.get("amplicons") + config.get("extra_regions")))

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


def compile_output_file_list(wildcards):
    outdir = pathlib.Path(output_spec.get("directory", "./"))
    output_files = []

    for f in output_spec["files"]:
        # Please remember to add any additional values down below
        # that the output strings should be formatted with.
        outputpaths = set(
            [
                f["output"].format(sample=sample, type=unit_type, target=target, report=report, caller=caller)
                for sample in get_samples(samples)
                for unit_type in get_unit_types(units, sample)
                for target in config.get("amplicons") + config.get("extra_regions")
                for report in config["multiqc"]["reports"]
                for caller in ["clairs_to", "deepsomatic", "vardict"]
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
        output_file = output_directory / pathlib.Path(f["output"])
        if pathlib.Path(f["output"]).suffix == '': # no file extension, hence a directory
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
                f'@workflow.output(directory("{output_file}"))' if pathlib.Path(f["output"]).suffix == '' else f'@workflow.output("{output_file}")',
                f'@workflow.log("logs/{rule_name}_{output_file.name}.log")',
                f'@workflow.container("{copy_container}")',
                f'@workflow.resources(time="{time}", threads={threads}, mem_mb="{mem_mb}", '
                f'mem_per_cpu={mem_per_cpu}, partition="{partition}")',
                f'@workflow.shellcmd("{copy_container}")', # replace with '@workflow.shellcmd("cp -r {input} {output}")\n\n' ?
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
