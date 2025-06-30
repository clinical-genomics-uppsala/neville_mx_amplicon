import pysam
import pandas as pd
import matplotlib.pyplot as plt

af_bins = [0.00,0.005,0.01,0.015,0.02,0.03,0.04,0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.5,0.6,0.7,0.8,0.9,1.0]

min_vaf = 0.1
vcf_file = f"/home/camille/ampliconthemato/pipeline_pool_amplicon/snv_indels/deepsomatic/min_vaf={min_vaf}/"\
           "D19-07233_T_reads.ont_adapt_trim.filtered.aligned.sorted.soft-clipped.deepsomatic.vcf.gz"
vcf = pysam.VariantFile(vcf_file)

varid = []
varvaf = []

for i, var in enumerate(vcf):
    varid.append(f"{var.chrom}:{var.pos}")
    varvaf.append(var.samples[0]['VAF'][0])

df_vaf = pd.DataFrame(zip(varid, varvaf), columns=['id', 'VAF']).set_index('id')

df_bins = pd.cut(df_vaf['VAF'], bins=af_bins, include_lowest=False)
df_bins.name = 'vaf_bin'
df_bincounts = df_bins.groupby(df_bins).count()
df_bincounts.name = 'vaf_bincounts'

fig, ax = plt.subplots(1, 1, figsize=(14, 10))
df_bincounts.plot(kind='bar', rot=90, ax=ax,
                  xlabel='VAF', ylabel='Count',
                  title=f'Counts of called variants per VAF-bin\n(total variants with VAF >= {min_vaf} = {df_vaf.shape[0]})')
plt.savefig(f"/home/camille/ampliconthemato/pipeline_pool_amplicon/snv_indels/deepsomatic/min_vaf={min_vaf}/histogram.png",
            bbox_inches='tight')
plt.show()

