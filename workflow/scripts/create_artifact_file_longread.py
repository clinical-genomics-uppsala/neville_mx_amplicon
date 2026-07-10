import gzip
import statistics
from collections import OrderedDict


vcf_files = snakemake.input.vcfs
artifact_panel_path = snakemake.output.artifact_panel
callers_used = snakemake.params.callers
print(f"Callers used by the pipeline: {callers_used}\n----------------------------------------------------\n")

# Dictionary with count of observations of a position and AF for each observation
call_stats_dict = {}


def get_info_field(info_str, key):
    """Safely extracts a value from the INFO field string."""
    for field in info_str.split(";"):
        if field.startswith(f"{key}="):
            return field.split("=", 1)[1]
    return None


def get_format_field(format_str, data_str, key):
    """Safely extracts a value from the FORMAT field string."""
    fields = data_str.split(":")
    keys = format_str.split(":")
    format_dict = dict(zip(keys, fields))
    return format_dict[key]


for file_name in vcf_files:
    print(f"call_stats_dict -- Processing {file_name}")
    with gzip.open(file_name, "rt") as infile:
        for line in infile:
            if line.startswith("#"):
                continue

            columns = line.strip().split("\t")
            if len(columns) < 8:
                continue

            chrom = columns[0]
            pos = int(columns[1])
            ref = columns[3]
            alt = columns[4]
            variant_type = "SNV" if (len(ref) == 1 and len(alt) == 1) else "INDEL"
            info = columns[7]
            formats = columns[8]
            data = columns[9]

            key = (chrom, pos, variant_type)

            callers_str = get_info_field(info, "CALLERS")
            print(f"\nVariant processed: {(chrom, pos, variant_type)} -- Callers {callers_str}")
            if not callers_str:
                continue
            callers = callers_str.split(",")
            print(f"\nVariant processed: {(chrom, pos, variant_type)} -- Callers SPLIT {callers}")

            af_str = get_format_field(formats, data, "AF")
            if not af_str:
                continue
            else:
                af = float(af_str)

            if key not in call_stats_dict:
                call_stats_dict[key] = {}  # {c: [0, []] for c in callers}

            for caller in callers:
                if caller not in call_stats_dict[key]:
                    call_stats_dict[key][caller] = [1, []]
                else:
                    call_stats_dict[key][caller][0] += 1
                    call_stats_dict[key][caller][1].append(af)

# Sort the variants per genomic position
call_stats_dict_sorted = OrderedDict(sorted(call_stats_dict.items()))
print('\n')
print(call_stats_dict_sorted)
print('\n')

with open(artifact_panel_path, "w") as artifact_panel:
    header = ["Chromosome", "pos", "variant_type"]
    for caller in callers_used:
        header.extend([caller, "median_MAF", "sd_MAF"])
    artifact_panel.write("\t".join(header) + "\n")

    for key, callers_data in call_stats_dict_sorted.items():
        chrom, pos, variant_type = key
        row = [str(chrom), str(pos), str(variant_type)]
        for caller in callers_used:
            obs_count = 0
            median_af = 0
            sd_af = 1000
            if caller not in callers_data:
                print(f"{caller} not in {callers_data}, extending row...\n")
                row.extend([f"{obs_count}", f"{median_af}", f"{sd_af}"])
            else:
                obs_count, af_list = callers_data[caller]
                af_list.sort()
                if len(af_list) >= 4:
                    print(f"Enough occurences of the position ({len(af_list)} observations), computing statistics for it...")
                    median_af = statistics.median(af_list)
                    sd_af = statistics.stdev(af_list)
                row.extend([f"{obs_count}", f"{median_af}", f"{sd_af}"])
        if sum([int(n) for n in row[3::3]]) == 0:
            print(f"\nNull counts for {row}")
            print(f"Data: {call_stats_dict[key]}\n")
        artifact_panel.write("\t".join(row) + "\n")
