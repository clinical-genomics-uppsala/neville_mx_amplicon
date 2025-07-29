__author__ = "Camille Clouard"
__copyright__ = "Copyright 2024, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"

import os

# Include pipeline specific rules
# include: "rules/common.smk"
# include: "rules/basecalling.smk"
# include: "rules/prefiltering.smk"
# include: "rules/aligning.smk"
# include: "rules/qc.smk"
# include: "rules/varcalling.smk"
# include: "rules/results_report.smk"

# 'All' rule, must be specified before any other modules are
# included, since they also contain 'All' rule
rule all:
    input:
        "annotation/bcbio_variation_recall_ensemble/D22-06934_T.ensembled.vep_annotated.view.vcf.gz"
        # unpack(compile_output_file_list),

# TODO: Draw DAG and rulegraph --> workflow/scripts
# snakemake -s workflow/Snakefile --configfile config/config.yaml --dag | dot -Tpdf > images/dag.pdf
# snakemake -s workflow/Snakefile --configfile config/config.yaml --rulegraph | dot -Tpdf > images/rulegraph.pdf

rule bcf_view_o_vcf:
    input:
        "annotation/bcbio_variation_recall_ensemble/D22-06934_T.ensembled.vep_annotated.vcf.gz",
    output:
        "annotation/bcbio_variation_recall_ensemble/D22-06934_T.ensembled.vep_annotated.view.vcf.gz",
    log:
        "annotation/bcbio_variation_recall_ensemble/D22-06934_T.ensembled.vep_annotated.view.vcf.gz.log",
    params:
        extra="-r chr17",
    container:
        "docker://hydragenetics/common:0.1.9"
    wrapper:
        "v1.24.0/bio/bcftools/view"