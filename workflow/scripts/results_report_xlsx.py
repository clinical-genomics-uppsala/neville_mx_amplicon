#!/bin/python3

from results_report_create_tables import *
from datetime import date
from operator import itemgetter
import subprocess
import yaml
import logging
import pandas as pd
import re

logging.basicConfig(
    format="{asctime} - {levelname} - {message}",
    style="{",
    datefmt="%Y-%m-%d %H:%M",
    level=logging.INFO,
)


def convert_columns_to_letter(nr_columns):
    # Function to convert number of columns to alphabetical coordinates for xlsx-sheets
    if nr_columns < 27:
        letter = chr(nr_columns + 64)
    elif nr_columns < 703:
        i = int((nr_columns - 1) / 26)
        letter = chr(i + 64) + chr(nr_columns - (i * 26) + 64)
    else:
        logging.error(f"Nr columns has to be less than 703, does not support three letter column-index for tables {nr_columns=}")
        sys.exit()
    return letter


vcf_snv_indels = snakemake.input.vcf_snv_indels
vcf_cnv_svs = snakemake.input.vcf_cnv_svs
sequenceid = snakemake.params.sequenceid
# poppy_version = snakemake.params.poppy_version

sample = snakemake.params.sample
sample_type = snakemake.params.sample_type
# Panels dict contain all info for sheets with subsets all SNVs, only need to change here if adding/removing subset sheets
# panels = {
#     "cll": {"bedfile": snakemake.input.cllbed, "vcf": snakemake.input.cll_vcf},
#     "myeloid": {"bedfile": snakemake.input.myeloidbed, "vcf": snakemake.input.myeloid_vcf},
# }

bedfiles = {
    "amplicons": snakemake.input.bedfile,
    "deepsomatic": snakemake.input.deepsomatic_bed,
    "vardict": snakemake.input.vardict_bed,
}

non_coding_regions = {
    "TERC": ["chr3", 169764300, 169766000, "entire gene + promotor"],
    "GATA2": ["chr3", 128481912, 128483849, "Intron4 (b/w exon4-5) in NM_001145662"],
    "ANKRD26": ["chr10", 27100074, 27100510, "Promotor and exon1"],
    "TP53": ["chr17", 7687366, 7687500, "Exon1 non-coding"],
    "NOTCH1": ["chr9", 136494400, 136496072, "3UTR"],
}

intron_coordinates = {}
for gene in non_coding_regions:
    chrom = non_coding_regions[gene][0]
    if chrom in intron_coordinates:  # If chrom in dict already
        intron_coordinates[chrom].append(non_coding_regions[gene][1:])
    else:
        intron_coordinates[chrom] = [non_coding_regions[gene][1:]]


synonymous_positions = {
    "NM_032638.5(GATA2):c.1416G>A": ["chr3", 128481046, "C", "T"],
    "NM_032638.5(GATA2):c.1023C>T": ["chr3", 128481939, "G", "A"],
    "NM_032638.5(GATA2):c.981G>A": ["chr3", 128483896, "C", "T"],
    "NM_032638.5(GATA2):c.649C>T": ["chr3", 128485949, "G", "A"],
    "NM_032638.5(GATA2):c.351C>G": ["chr3", 128486247, "G", "C"],
    "NM_000546.6(TP53):c.375G>A": ["chr17", 7675994, "C", "T"],
    "NM_000546.6(TP53):c.375G>T": ["chr17", 7675994, "C", "A"],
    "NM_000546.6(TP53):c.375G>C": ["chr17", 7675994, "C", "G"],
    "NM_000546.6(TP53):c.672G>A": ["chr17", 7674859, "C", "T"],
    "NM_000546.6(TP53):c.993G>A": ["chr17", 7673535, "C", "T"],
}


# wanted_transcripts = []
# with open(snakemake.input.wanted_transcripts) as wanted_file:
#     for line in wanted_file:
#         wanted_transcripts.append(line.split()[1].split(".")[0])

# thresholds = [int(x) for x in snakemake.params.thresholds.split(",")]

for x in VariantFile(vcf_snv_indels).header.records:
    if x.key == "VEP":
        vep_line = x.value

# Add filters from yaml-files
filters_snv = []
filters_svs = []

for filter_file in [snakemake.params.filter_somatic_hard, snakemake.params.filter_somatic_soft]:
    with open(filter_file, "r") as file:
        filters_dict = yaml.safe_load(file)

    for key, items in filters_dict["filters"].items():
        if "hard_filter_somatic" in filter_file:
            filters_snv.append(key + ": " + items["description"])
            # filters_p.append(key + ": " + items["description"])
        elif "soft_filter_somatic" in filter_file:
            filters_snv.append(items["soft_filter_flag"] + ": " + items["description"])
        # elif "soft_filter_p" in filter_file:
        #     filters_p.append(items["soft_filter_flag"] + ": " + items["description"])
        else:
            logging.error(f"Unknown filter file, have the filename changed? {filter_file=}")


""" Create data tables to populate excel sheets with """
logging.info(f"Creating data tables to populate excel sheets with")
var_table = VariantTable(vcf_snv_indels=vcf_snv_indels, vcf_cnv_sv=vcf_cnv_svs)
var_table.create_snv_table(sequenceid=sequenceid)
var_table.create_sv_table(sequenceid=sequenceid)

snv_table = var_table.get_table_snv_indels()
structv_table = var_table.get_table_cnv_sv()
# known_table, known_percent = create_known_variants_table(vcf, vcf_cnv_svs, sequenceid)

# for panel in panels.keys():
#     logging.debug(f"Creating {panel} table")
#     panels[panel]["table"] = create_snv_table(panels[panel]["vcf"], sequenceid)
# logging.debug(f"{panels=}")


logging.debug(f"Intron and synonymous tables")
intron_table = []
synonymous_table = []
for record in snv_table["data"]:
    if record[4] in intron_coordinates:
        for pair in intron_coordinates[record[4]]:
            if record[5] >= pair[0] and record[5] <= pair[1]:
                intron_table.append(record)

    for position_list in [coordinate for coordinate in synonymous_positions.values() if record[4] == coordinate[0]]:
        if position_list[1] == record[5] and position_list[2] == record[6] and position_list[3] == record[7]:
            synonymous_table.append(record)  # AF filter?


# logging.debug(f"List genes in  bed")
# _genes = []
# with open(snakemake.input.bed, "r") as _file:
#     for lline in _file:
#         line = lline.strip().split("\t")
#         _genes.append(line[3])
# _genes = list(dict.fromkeys(_genes))

logging.debug(f"Avg cov per amplicon (including overlapping regions) from mosdepth {snakemake.input.mosdepth_regions=}")
regionscov_table = {"data": [], "headers": []}
regionscov_table["headers"] = [
    {"header": "Chr"},
    {"header": "Start"},
    {"header": "Stop"},
    {"header": "Gene"},
    {"header": "Target"},
    {"header": "Avg Coverage"},
]
bed_table = []
with (gzip.open(snakemake.input.mosdepth_regions, "rt") as regions_file):
    for lline in regions_file:
        line = lline.strip().split("\t")
        gene = line[3].split("_")[0]
        # if gene == "FLT3":
        #     gene += "-" + line[3].split("_")[1]
        # transcript = ""  # "_".join(line[3].split("_")[1:3])
        # exon = ""  # str(line[3].split("_")[3])
        target = line[3].replace("_3kb", "").split("_")[1]
        if re.search(r"^[A-Z]+.*", target) is None:
            target = ""
        coverage_row = [line[0], line[1], line[2], gene, target, float(line[4])]
        if coverage_row not in regionscov_table["data"]:
            regionscov_table["data"].append(coverage_row)
        if line[0:5] not in bed_table:
            bed_table.append(line[0:5])
# TODO: Add caller-specific BED info

logging.debug(f"Estimated reads counts per amplicon (avg coverage excluding overlapping regions) from mosdepth {snakemake.input.mosdepth_regions=}")
targetcounts_table = {"data": [], "headers": []}
targetcounts_table["headers"] = [
    {"header": "Amplicon"},  # , "total_string": "Total counts"},
    {"header": "Estimated Reads Counts"}  # , "total_function": "sum"},
]
targets = [a for a in snakemake.params.amplicons if a.find("TP53") < 0] + snakemake.params.extra_regions  # select the non overlapping regions
csvcounts = pd.read_csv(snakemake.input.csv_counts, index_col="target")  #, sep="\t")
csvcounts.index = csvcounts.index.str.replace(f"{sample}_{sample_type}_", "")
tot_reads = 0
for counts_row in csvcounts.itertuples():
    if counts_row[0] in targets:
        tot_reads += round(float(counts_row[-1]))
        if counts_row[0].find("D2+J3") < 0:
            targetcounts_table["data"].append([counts_row[0].replace("_only", ""),
                                               round(float(counts_row[-1]))
                                               ]
                                              )
        else:  # J3 = (D2+J3) - D2
            countsD = csvcounts.loc["TP53_D2_only"]["mean"]
            targetcounts_table["data"].append(["TP53_J3",
                                               round(counts_row[-1] - countsD)
                                               ]
                                              )
targetcounts_table["data"].append(["Total counts", tot_reads])

logging.debug(f"Estimated reads counts per amplicon in each pool (avg coverage excluding overlapping regions)")
poolcounts_table = []
for p, csv_pool in enumerate(snakemake.input.pool_counts):
    poolcounts = pd.read_csv(csv_pool)
    poolcounts.set_index(poolcounts.columns[0], inplace=True)
    poolcounts_table.append({"data": [], "headers": []})
    poolcounts_table[-1]["headers"] = [
        {"header": poolcounts.index.name},
        {"header": "Estimated Reads Counts"},
        {"header": "% of Reads Counts"},
    ]
    poolcounts_table[-1]["data"] = [[counts_row[0].replace(f"{sample}_{sample_type}_", "").replace("_only", ""),
                                     *counts_row[1:]] for counts_row in poolcounts.itertuples()]
    poolcounts_table[-1]["data"].append(["", ""])  # leave blank line
    poolcounts_table[-1]["data"].append(["% of Seq. Run", round(100 * poolcounts.loc["total"].iloc[0] / tot_reads, 1)])


logging.debug(f"Overview sheets qc-table values extracted")
coverage = {}
with open(snakemake.input.mosdepth_summary, "r") as summary_file:
    for lline in summary_file:
        line = lline.strip().split("\t")
        if line[0] == "total_region":
            coverage["avg_cov"] = line[3]
        # elif line[0] == "chrX_region":
        #     coverage["chrX_cov"] = line[3]
        # elif line[0] == "chrY_region":
        #     coverage["chrY_cov"] = line[3]

logging.debug(f"Avg cov per exon/coding region from mosdepth {snakemake.input.mosdepth_exons=}")
exonscov_table = {"data": [], "headers": []}
exonscov_table["headers"] = [
    {"header": "Chr"},
    {"header": "Start"},
    {"header": "Stop"},
    {"header": "Gene"},
    {"header": "Exon"},
    {"header": "Transcript"},
    {"header": "Avg Coverage"},
]
bed_table = []
with gzip.open(snakemake.input.mosdepth_exons, "rt") as exons_file:
    for lline in exons_file:
        line = lline.strip().split("\t")
        print(line)
        gene = line[3].split("_")[0]
        transcript = "_".join(line[3].split("_")[1:3])
        exon = "" # str(line[3].split("_")[3])
        coverage_row = [line[0], line[1], line[2], gene, exon, transcript, float(line[4])]
        if coverage_row not in regionscov_table["data"]:
            exonscov_table["data"].append(coverage_row)
        if line[0:5] not in bed_table:
            bed_table.append(line[0:5])


""" xlsx file with sheets """
logging.info(f"Creating xlsx file {snakemake.output.xlsx}")
workbook = xlsxwriter.Workbook(snakemake.output.xlsx)

worksheet_overview = workbook.add_worksheet("Overview")
# if sample.lower() == "hd829":
#     worksheet_known = workbook.add_worksheet("Known variants")
# else:
#     for panel in panels.keys():
#         panels[panel]["sheet"] = workbook.add_worksheet(panel.upper())
worksheet_snv = workbook.add_worksheet("SNV_Indels")
worksheet_sv = workbook.add_worksheet("SVs")
# worksheet_intron = workbook.add_worksheet("Intron")
# worksheet_syno = workbook.add_worksheet("Synonymous")
worksheet_cov = workbook.add_worksheet("Coverage")
worksheet_exonscov = workbook.add_worksheet("Exon Coverage")


empty_list = ["", "", "", "", "", ""]
format_heading = workbook.add_format({"bold": True, "font_size": 18})
format_line = workbook.add_format({"top": 1})
format_bold = workbook.add_format({"bold": True})
format_orange = workbook.add_format({"bg_color": "#ffd280"})
format_red = workbook.add_format({"font_color": "red"})
format_table_heading = workbook.add_format({"bold": True, "text_wrap": True})
float_format = workbook.add_format({"num_format": "# ##0.00"})

""" Overview """
logging.debug(f"Overview sheet")
worksheet_overview.write(0, 0, sample, format_heading)
worksheet_overview.write(1, 0, "RunID: " + sequenceid)
worksheet_overview.write(2, 0, "Processing date: " + date.today().strftime("%B %d, %Y"))
worksheet_overview.write_row(3, 0, empty_list, format_line)

worksheet_overview.write(4, 0, "Created by: ")
worksheet_overview.write(4, 4, "Valid from: ")
worksheet_overview.write(5, 0, "Signed by: ")
worksheet_overview.write(5, 4, "Document nr: ")
worksheet_overview.write_row(6, 0, empty_list, format_line)

worksheet_overview.write(7, 0, "Sheets:", format_table_heading)
i = 8
# if sample.lower() == "hd829":
#     worksheet_overview.write_url(i, 0, "internal:'Known variants'!A1", string="Known variants")
# else:
#     for panel in panels.keys():
#         worksheet_overview.write_url(i, 0, "internal: '" + panel.upper() + "'!A1", string=panel.upper() + " panel variants")
#         i += 1
worksheet_overview.write_url(i, 0, "internal:'SNV_Indels'!A1", string="SNVs and Indels identified")
worksheet_overview.write_url(i + 1, 0, "internal:'SVs'!A1", string="SVs identified")
# worksheet_overview.write_url(i + 2, 0, "internal: 'Intron'!A1", string="Intron and non-coding variants")
# worksheet_overview.write_url(i + 3, 0, "internal: 'Synonymous'!A1", string="Synonymous variants")
# worksheet_overview.write_url(i + 4, 0, "internal:'Low Coverage'!A1", string="Low Coverage regions")
worksheet_overview.write_url(i + 5, 0, "internal:'Coverage'!A1", string="Coverage")
# worksheet_overview.write_url(i + 6, 0, "internal: 'QCI'!A1", string="QCI")
i += 8

worksheet_overview.write_row(
    i,
    0,
    [
        "RunID",
        "DNAnr",
        "Avg. coverage [x]",
        # "Duplicationlevel [%]",
        # str(thresholds[0]) + "x",
        # str(thresholds[1]) + "x",
        # str(thresholds[2]) + "x",
    ],
    format_table_heading,
)
worksheet_overview.write_row(i + 1, 0, [sequenceid, sample, coverage["avg_cov"]])  # , str(duplication_rate)] + thresholds_results)
i += 3

worksheet_overview.write(i, 0, "Average coverage of regions in amplicons bedfile")
# worksheet_overview.write_row(i + 1, 0, ["chrX", coverage["chrX_cov"]])
# worksheet_overview.write_row(i + 2, 0, ["chrY", coverage["chrY_cov"]])
i += 4

# worksheet_overview.write(i, 0, "Twist Myeloid capture panel TE-96463385_hg38 used")
worksheet_overview.write(
    i + 1,
    0,
    "with the pipeline pipeline_pool_amplicon.",
)
worksheet_overview.write_url(i + 3, 0, "https://github.com/camcl/pipeline_pool_amplicon", string="Pipeline repository")

worksheet_overview.write(i + 6, 0, "Specific program versions can be found in MultiQC report")
i += 8

worksheet_overview.write(i, 0, "Full design bedfile: " + bedfiles["amplicons"])
worksheet_overview.write(i + 1, 0, "DeepSomatic-specific bedfile (exons in TP53, IDH1, and IDH2): " + bedfiles["deepsomatic"])
worksheet_overview.write(i + 2, 0, "VarDict-specific bedfile (ITDs in FLT3 and NPM1): " + bedfiles["vardict"])
# worksheet_overview.write(i + 3, 0, "Background panel used for snvs: " + snakemake.params.background)
# worksheet_overview.write(i + 4, 0, "indel bedfile used: " + snakemake.input.indelbed)
# worksheet_overview.write(i + 5, 0, "indel artifact panel used: " + snakemake.params.artifact_indel)
i += 6


""" SNVs sheet """
logging.debug(f"SNV_Indels sheet")
worksheet_snv.set_column(2, 2, 10)
worksheet_snv.set_column(5, 5, 10)
worksheet_snv.set_column(11, 13, 10)
worksheet_snv.write("A1", "Variants found", format_heading)
worksheet_snv.write("A3", "Sample: " + str(sample))
worksheet_snv.write("A4", "Reference used: " + str(snakemake.params.ref))
worksheet_snv.write("A6", "Databases used: " + vep_line)

worksheet_snv.write("A8", "Filters: ", format_orange)
for i, filter_txt in enumerate(filters_snv):
    i += 9
    worksheet_snv.write("B" + str(i), filter_txt, format_orange)

i += 2
worksheet_snv.write_rich_string(
    "A" + str(i), "Only variants with ", format_bold, "> 3 % AF", " and filter-flag ", format_bold, "PASS", " shown by default."
)
worksheet_snv.write(
    "A" + str(i + 1),
    "To see all variants; put marker on header row, then click on 'Standard Filter' and remove any values. "
    + "You can then use the drop-downs in the header row to filter to your liking.",
)

# worksheet_snv.write_url(
#     "A" + str(i + 3),
#     "external:" + sample + "_" + sample_type + "_" + sequenceid + "_bamsnap/index.html",
#     string="SNV screenshots",
# )
# worksheet_snv.write("A" + str(i + 4), "Only variants with AF >= 5% and PASS have automated screenshots.")

i += 6
column_end = ":" + convert_columns_to_letter(len(snv_table["headers"]))
if len(snv_table["data"]) > 0:
    table_area = "A" + str(i) + column_end + str(len(snv_table["data"]) + i)
    table_area_data = "A" + str(i + 1) + column_end + str(len(snv_table["data"]) + i)
else:
    table_area = "A" + str(i) + column_end + str(i + 1)
    table_area_data = "A" + str(i + 1) + column_end + str(i + 1)

worksheet_snv.add_table(table_area, {"columns": snv_table["headers"],
                                     "style": "Table Style Light 1",
                                     "autofilter": False
                                     })
cond_formula = "=LEFT($A" + str(i + 1) + ', 4)<>"PASS"'
worksheet_snv.conditional_format(table_area_data, {"type": "formula", "criteria": cond_formula, "format": format_orange})


worksheet_snv.autofilter(table_area)
worksheet_snv.filter_column("A", "Filter != PASS")
worksheet_snv.filter_column("I", "AF >= 0.03")
for row_data in snv_table["data"]:
    if row_data[0] == "PASS" and float(row_data[8]) >= 0.03:
        pass
    else:
        worksheet_snv.set_row(i, options={"hidden": True})
    worksheet_snv.write_row(i, 0, row_data)
    i += 1

if True:
    """ SVs sheet """
    logging.debug(f"sv sheet")
    worksheet_sv.set_column(2, 2, 10)
    worksheet_sv.set_column(5, 5, 10)
    worksheet_sv.set_column(11, 13, 10)
    worksheet_sv.write("A1", "Variants found", format_heading)
    worksheet_sv.write("A3", "Sample: " + str(sample))
    worksheet_sv.write("A4", "Reference used: " + str(snakemake.params.ref))
    # worksheet_sv.write("A6", "To limit runtime indel were used with a specific designfile: " + snakemake.input.indelbed)
    # worksheet_sv.write("A7", "Which includes the following genes: ")
    i = 8
    # for gene in indel_genes:
    #     worksheet_sv.write("C" + str(i), gene)
    #     i += 1

    worksheet_sv.write("A" + str(i + 1), "Filters: ", format_orange)
    # filters_svs = filters_snv
    for j, filter_txt in enumerate(filters_svs):
        j += i + 1
        worksheet_sv.write("B" + str(j), filter_txt, format_orange)
    i += 2 + len(filters_svs)

    worksheet_sv.write_rich_string("A" + str(i), "Only variants with filter-flag ", format_bold, "PASS", " shown by default.")
    worksheet_sv.write(
        "A" + str(i + 1),
        "To see all variants; put marker on header row, then click on 'Standard Filter' and remove any values. "
        + "You can then use the drop-downs in the header row to filter to your liking.",
    )

    i += 3
    column_end = ":" + convert_columns_to_letter(len(structv_table["headers"]))
    if len(structv_table["data"]) > 0:
        structv_table_area = "A" + str(i) + column_end + str(len(structv_table["data"]) + i)
        structv_table_area_data = "A" + str(i + 1) + column_end + str(len(structv_table["data"]) + i)
    else:
        structv_table_area = "A" + str(i) + column_end + str(i + 1)
        structv_table_area_data = "A" + str(i + 1) + column_end + str(i + 1)

    worksheet_sv.add_table(structv_table_area,
                           {"columns": structv_table["headers"],
                            "style": "Table Style Light 1",
                            "autofilter": False
                            })
    cond_formula = "=LEFT($A" + str(i + 1) + ', 4)<>"PASS"'
    worksheet_sv.conditional_format(structv_table_area_data, {"type": "formula", "criteria": cond_formula, "format": format_orange})

    worksheet_sv.autofilter(structv_table_area)
    worksheet_sv.filter_column("A", "Filter != PASS")
    for row_data in structv_table["data"]:
        if row_data[0] == "PASS":
            pass
        else:
            worksheet_sv.set_row(i, options={"hidden": True})
        worksheet_sv.write_row(i, 0, row_data)
        i += 1


""" Coverage sheet"""
logging.debug(f"Coverage sheet")
worksheet_cov.set_column(1, 2, 10)
worksheet_cov.set_column(5, 5, 15)
worksheet_cov.write(0, 0, "Average Coverage per amplicon (including overlapping regions)", format_heading)
worksheet_cov.write_row(1, 0, empty_list, format_line)
worksheet_cov.write(2, 0, "Sample: " + str(sample))
worksheet_cov.write(3, 0, "Average coverage of each region in amplicons-bedfile")

top_row = 6
tab_margin = 3  # leave 3-1=2 blank lines/columns between tables

"""Table coverage per amplicon"""
cov_tab_params = {"first_col": "A",
                  "last_col": convert_columns_to_letter(len(regionscov_table["headers"])),
                  "width": len(regionscov_table["headers"]),
                  "height": len(regionscov_table["data"]),
                  "offset_top": top_row
                  }

cov_table_area = f"{cov_tab_params['first_col']}{cov_tab_params['offset_top']}" \
                 + ":" \
                 + f"{cov_tab_params['last_col']}{cov_tab_params['height'] + cov_tab_params['offset_top']}"

worksheet_cov.add_table(
    cov_table_area, {"data": regionscov_table["data"],
                     "columns": [{"header": hd["header"], "format": float_format} for hd in regionscov_table["headers"]],
                     "style": "Table Style Light 1",
                     "autofilter": False
                     }
)

"""Table estimated counts per amplicon"""
counts_tab_params = {"first_col": "A",
                     "last_col": convert_columns_to_letter(len(targetcounts_table["headers"])),
                     "width": len(targetcounts_table["headers"]),
                     "height": len(targetcounts_table["data"]),
                     "offset_top": top_row + cov_tab_params['height'] + tab_margin  # place the table below coverage
                     }

counts_table_area = f"{counts_tab_params['first_col']}{counts_tab_params['offset_top']}" \
    + ":" \
    + f"{counts_tab_params['last_col']}{counts_tab_params['height'] + counts_tab_params['offset_top']}"

worksheet_cov.add_table(
    counts_table_area, {"data": targetcounts_table["data"],
                        "columns": targetcounts_table["headers"],
                        "style": "Table Style Light 1",
                        "autofilter": False
                        }
                        # "total_row": True}
)

"""Tables reads counts per amplicon in each pool"""
for p, tab_pool in enumerate(poolcounts_table):
    print(tab_pool)
    pools_tab_params = {"first_col": convert_columns_to_letter(1 + p * (tab_margin + len(tab_pool["headers"]))),  # letter A = column 1 and letter @ = column 0
                        "last_col": convert_columns_to_letter(p * (tab_margin + len(tab_pool["headers"])) + len(tab_pool["headers"])),
                        "width": len(tab_pool["headers"]),
                        "height": len(tab_pool["data"]),
                        "offset_top": top_row + cov_tab_params['height'] + counts_tab_params['height'] + 2 * tab_margin  # place the table below coverage and estimated reads counts
                        }

    pools_table_area = f"{pools_tab_params['first_col']}{pools_tab_params['offset_top']}" \
        + ":" \
        + f"{pools_tab_params['last_col']}{pools_tab_params['height'] + pools_tab_params['offset_top']}"
    print(pools_table_area)

    worksheet_cov.add_table(
        pools_table_area, {"data": tab_pool["data"],
                           "columns": tab_pool["headers"],
                           "style": "Table Style Light 1",
                           "autofilter": False
                           }
    )

"""Lineplot with sequencing throughput"""
worksheet_cov.insert_image(f"{convert_columns_to_letter(len(regionscov_table['headers']) + tab_margin)}{top_row}",
                           snakemake.input.yield_plot,
                           {"x_scale": 0.65, "y_scale": 0.85}
                           )

""" Exon coverage sheet"""
logging.debug(f"Exon coverage sheet")
worksheet_exonscov.set_column(1, 2, 10)
worksheet_exonscov.set_column(5, 5, 15)
worksheet_exonscov.write(0, 0, "Average Coverage per Exon", format_heading)
worksheet_exonscov.write_row(1, 0, empty_list, format_line)
worksheet_exonscov.write(2, 0, "Sample: " + str(sample))
worksheet_exonscov.write(3, 0, "Average coverage of each region in exon-bedfile")

column_end = ":" + convert_columns_to_letter(len(exonscov_table["headers"]))
table_area = "A6" + column_end + str(len(exonscov_table["data"]) + 6)

worksheet_exonscov.add_table(
    table_area, {"data": exonscov_table["data"],
                 "columns": exonscov_table["headers"],
                 "style": "Table Style Light 1",
                 "autofilter": False
                 }
)

workbook.close()
logging.info(f"All done!")
