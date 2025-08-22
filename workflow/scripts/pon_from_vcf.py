# Adapted from https://github.com/KolmogorovLab/Severus/blob/main/scripts/pon_from_vcf.py

# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  9 14:25:41 2024

@author: keskusa2
"""

import pysam

vcf_file = '/home/camille/ampliconthemato/pipeline_pool_amplicon/reference/pon_CU11_C267.vcf.gz'
vcf = pysam.VariantFile(vcf_file)
n_samples = 2
sv_list = []
for var in vcf:
    print(*var.info)
    cc = sum([1 for idd in list(var.samples) if not var.samples[idd]['GT'] == (None, None)]) / n_samples
    print(cc)
    print(var.chrom, str(var.pos), var.chrom)
    if var.info['SVTYPE'] == 'INS':
        sv_list.append(','.join([var.chrom, str(var.pos), var.chrom, str(var.info['SVLEN']), str(var.info['CIPOS'][1]),
                                 str(var.info['CIEND'][1]), var.info['SVTYPE'], str(cc)]))
    elif var.info['SVTYPE'] == 'BND':
        sv_list.append(','.join([var.chrom, str(var.pos), var.info['CHR2'], str(var.stop), str(var.info['CIPOS'][1]),
                                 str(var.info['CIEND'][1]), var.info['SVTYPE'], str(cc)]))
    else:
        sv_list.append(','.join(
            [var.chrom, str(var.pos), var.chrom, str(var.stop), str(var.info['CIPOS'][1]), str(var.info['CIEND'][1]),
             var.info['SVTYPE'], str(cc)]))

out_file3 = '/home/camille/ampliconthemato/pipeline_pool_amplicon/reference/PON.tsv'
with open(out_file3, "w") as fout3:
    for line in sv_list:
        fout3.write(line)
        fout3.write('\n')
fout3.close()