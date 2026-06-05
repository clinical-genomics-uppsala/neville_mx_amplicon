import pandas as pd
import pathlib
from pathlib import Path
from snakemake.utils import validate
from snakemake.utils import min_version
import sys
import yaml

from hydra_genetics.utils.resources import load_resources
from hydra_genetics import min_version as hydra_min_version
from hydra_genetics.utils.misc import replace_dict_variables
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *


hydra_min_version("3.0.0")
min_version("7.32.0")

if not workflow.overwrite_configfiles:
    sys.exit("config file has to be specified with --configfile")

config = replace_dict_variables(config)

# Validate config
# validate(config, schema="../schemas/config_references.schema.yaml")
config = load_resources(config, config["resources"])
config = load_resources(config, config["resources_references"])
validate(config, schema="../schemas/resources.schema.yaml")

# Sample information
samples = pd.read_table(config["samples_ref"], comment="#").set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units_ref"], dtype=str).set_index(["sample", "type"], drop=False)
panel_ref =  pd.read_table(config["panel_ref"], dtype=str).set_index(["experiment", "sample", "type"], drop=False)
validate(units, schema="../schemas/units.schema.yaml")

# Ouput file specification
with open(config["output"]) as f:
    output_spec = yaml.safe_load(f)
    validate(output_spec, schema="../schemas/output_files.schema.yaml")

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


def get_gvcfs_clairs_to():
    return list(set([f"reference_files/{t.experiment}_{t.sample}_{t.type}_snv.clairs_to.g.vcf.gz"
                     for t in panel_ref.itertuples()
                     ])
    )


def get_gvcfs_deepsomatic():
    return list(set([f"reference_files/{t.experiment}_{t.sample}_{t.type}.deepsomatic.g.vcf.gz"
                     for t in panel_ref.itertuples()
                     ])
    )


def get_gvcfs_vardict():
    return list(set([f"reference_files/{t.experiment}_{t.sample}_{t.type}.vardict.g.vcf.gz"
                     for t in panel_ref.itertuples()
                     ])
    )


def get_vcfs():
    return list(
        set(
            [
                f"snv_indels/bcbio_variation_recall_ensemble/{t.experiment}_{t.sample}_{t.type}.ensembled.vep_annotated.rename_vaf.filter.somatic_soft.vcf.gz"
                for t in panel_ref.itertuples()
            ]
        )
    )