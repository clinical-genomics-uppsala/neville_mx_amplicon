# vim: syntax=python tabstop=4 expandtab
# coding: utf-8

__author__ = "Camille Clouard"
__copyright__ = "Copyright 2024, Camille Clouard"
__email__ = "camille.clouard@scilifelab.uu.se"
__license__ = "GPL-3"


import pysam
def dummy():
    return 1

def get_annotation_data_format(sample_index):
    return lambda variant, field: variant.samples[sample_index].get(field, None)

vcf = pysam.VariantFile("/home/camille/ampliconthemato/pipeline_pool_amplicon/annotation/bcbio_variation_recall_ensemble/D22-06934_T.ensembled.vep_annotated.vcf.gz")

for rec in vcf:
    print(*zip(rec.samples['D22-06934'].keys()), rec.samples['D22-06934'].get('AF'))
    print(get_annotation_data_format('D22-06934'))