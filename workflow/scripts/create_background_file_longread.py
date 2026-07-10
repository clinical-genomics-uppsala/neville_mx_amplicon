import gzip
import statistics
from collections import OrderedDict


gvcf_files = snakemake.input.gvcfs
background_file_path = snakemake.output.background_file
min_dp = snakemake.params.min_dp
max_af = snakemake.params.max_af

background_dict = {}

for file_name in gvcf_files:
    with gzip.open(file_name, "rt") as infile:
        for line in infile:
            if line.startswith("#"):
                continue

            columns = line.strip().split("\t")
            if len(columns) < 10:
                continue

            # Handle potential chromosome prefixes
            chrom = columns[0]
            if chrom.startswith("chr"):
                chrom = chrom[3:]

            pos = columns[1]
            key = (chrom, pos)

            formats = columns[8].split(":")
            data = columns[9].split(":")

            try:
                ad_idx = formats.index("AD")
                ad_info = data[ad_idx].split(",")
                # handle specific format in ClairS-TO for RefCalls, for example:
                # GT:GQ:DP:AF:AD:AU:CU:GU:TU   0/0:8:7181:0.9109:6541:2:6541:1:0
                # --> count for alt alleles are :-separated
                if len(ad_info) == 1:
                    ad_info = data[ad_idx:]  # AD:AU:CU:GU:TU    6541:2:6541:1:0
                    print(ad_info)

                ref_ad = int(ad_info[0])  # 6541
                alt_ad = sum(int(ad) for ad in ad_info[1:])  # sum([2, 6541, 1, 0]) = 6544
                # or alt_ad = sum(int(ad) for ad in ad_info[1:]) - int(max(ad_info[1:]))  to sum only the alternate alleles
                dp = ref_ad + alt_ad

                if dp <= min_dp:
                    continue

                alt_af = alt_ad / float(dp)  # 1.0005

                # Correct VAF for RefCalls in ClairS-TO
                if alt_af > 1:
                    alt_af = alt_af - 1  # 0.0005 < max_af=0.03

                # Check if it's potentially a homozygous variant being filtered
                if alt_af > 1 - max_af:  # 1.0005 > 1 - 0.03
                    alt_af = 1 - alt_af

                if alt_af > max_af:
                    continue

                if key in background_dict:
                    background_dict[key].append(alt_af)
                else:
                    background_dict[key] = [alt_af]
            except (ValueError, IndexError):
                # Skip if AD field or index is missing/malformed
                continue

# Sort the variants per genomic position
background_dict_sorted = OrderedDict(sorted(background_dict.items()))

# Write background_dict to json file
# json_file_path = snakemake.output.background_json
# import json
# with open(json_file_path, 'w') as json_file:
#     json.dump(background_dict_sorted, json_file)

with open(background_file_path, "w") as background_file:
    background_file.write("Chr\tPos\tMedian_or_Mean\tSD\tNrObs\n")
    for key, af_list in background_dict_sorted.items():
        if len(af_list) < 4 and len(af_list) > 0:
            median_background = statistics.mean(af_list)
            stdev_background = 0
        else:
            median_background = statistics.median(af_list)
            stdev_background = statistics.stdev(af_list)

        chrom, pos = key
        background_file.write(f"{chrom}\t{pos}\t{median_background}\t{stdev_background}\t{len(af_list)}\n")
