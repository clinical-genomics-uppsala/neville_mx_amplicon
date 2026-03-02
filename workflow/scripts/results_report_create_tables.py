#!/bin/python3

import xlsxwriter
import gzip
from pysam import VariantFile


# VEP fields in list to get index
def index_vep(variantfile):
    csq_index = []
    for x in variantfile.header.records:
        if "CSQ" in str(x):
            csq_index = str(x).split("Format: ")[1].strip().strip('">').split("|")
    return csq_index


# Extract table columns from vcf records
def extract_vcf_values(record, csq_index):
    return_dict = {}
    csq = record.info["CSQ"][0].split("|")
    sample = record.samples.keys()[0]
    variant = record.samples.values()[0]

    try:
        return_dict["callers"] = ";".join(record.info["CALLERS"])
        caller = record.info["CALLERS"][0]
    except KeyError:
        return_dict["callers"] = "Sniffles2"
        caller = "Sniffles2"

    print("\n-------------------CALLER:", caller)

    if caller == "vardict":
        try:
            return_dict["af"] = float(record.info["AF"][0])
        except KeyError:
            return_dict["af"] = int(record.samples.values()[0].get("AD")[1]) / sum(record.samples.values()[0].get("AD"))

        return_dict["dp"] = record.samples.values()[0]["DP"]

    elif caller in ["clairs_to", "deepsomatic"]:
        try:
            return_dict["af"] = float(record.samples.values()[0].get("AF")[0])
        except KeyError:
            return_dict["af"] = 999.0

        return_dict["dp"] = record.samples.values()[0]["DP"]

    else:  # Sniffles2
        try:
            return_dict["af"] = float(record.info["VAF"])
            return_dict["af"] = int(record.info["SUPPORT"]) / sum([int(cov) for cov in record.info["COVERAGE"]])
        except KeyError:
            return_dict["af"] = 999.0

        return_dict["dp"] = sum([int(cov) for cov in record.info["COVERAGE"]]) / len(record.info["COVERAGE"])

    try:
        return_dict["svlen"] = int(record.info["SVLEN"])
    except KeyError:
        pass

    try:
        return_dict["artifact_callers"] = (
            str(record.info["Artifact"])
            .replace("(", "")
            .replace(")", "")
            .replace(", ", ";")
            .replace("'", "")
        )
        if type(record.info["ArtifactMedian"]) == str:
            return_dict["artifact_median"] = str(round(float(record.info["ArtifactMedian"]), 3))
        else:
            return_dict["artifact_median"] = ";".join([str(round(float(x), 3)) for x in record.info["ArtifactMedian"]])
        return_dict["artifact_nr_sd"] = (str(record.info["ArtifactNrSD"])
                                         .replace("(", "")
                                         .replace(")", "")
                                         .replace(", ", ";"))
    except KeyError:
        pass

    try:
        return_dict["background_median"] = record.info["PanelMedian"]
        return_dict["background_nr_sd"] = record.info["PositionNrSD"]
    except KeyError:
        return_dict["background_median"] = ""
        return_dict["background_nr_sd"] = ""

    return_dict["gene"] = csq[csq_index.index("SYMBOL")]
    return_dict["transcript"] = csq[csq_index.index("HGVSc")].split(":")[0]

    if len(csq[csq_index.index("HGVSc")].split(":")) > 1:
        return_dict["coding_name"] = csq[csq_index.index("HGVSc")].split(":")[1]
    else:
        return_dict["coding_name"] = ""
    return_dict["ensp"] = csq[csq_index.index("HGVSp")]
    return_dict["consequence"] = csq[csq_index.index("Consequence")]

    existing = csq[csq_index.index("Existing_variation")].split("&")
    cosmic_list = [cosmic for cosmic in existing if cosmic.startswith("CO")]
    if len(cosmic_list) == 0:
        return_dict["cosmic"] = ""
    else:
        return_dict["cosmic"] = ", ".join(cosmic_list)

    return_dict["clinical"] = csq[csq_index.index("CLIN_SIG")]

    rs_list = [rs for rs in existing if rs.startswith("rs")]
    if len(rs_list) == 0:
        return_dict["rs"] = ""
    else:
        return_dict["rs"] = ", ".join(rs_list)
    return_dict["max_pop_af"] = csq[csq_index.index("MAX_AF")]
    return_dict["max_pops"] = csq[csq_index.index("MAX_AF_POPS")]
    return_dict["filter_flag"] = ",".join(record.filter.keys())
    if return_dict["af"] >= 0.05:
        return_dict["igv"] = "pathpath"
    else:
        return_dict["igv"] = ""

    return return_dict


# def create_known_variants_table(vcf_input, pindel_input, sequenceid):
#     known_number = 0
#     vcf_file = VariantFile(vcf_input)
#     pindel_file = VariantFile(pindel_input)
#     sample = list(vcf_file.header.samples)[0]
#     known = [
#         ["IDH1", "R132C", "c.394C>T", "chr2", 208248389, "G", "A", "COSM28747", "COSV61615256", 0.0442],
#         ["IDH2", "R172K", "c.515G>A", "chr15", 90088606, "C", "T", "COSM33733", "COSV57468734", 0.0262],
#         ["NPM1", "W288Cfs*12", "c.860_863dup", "chr5", 171410539, "C", "CTCTG", "COSM17559", "COSV51542664", 0.0389],
#         ["FLT3", "D835Y", "c.2503G>T", "chr13", 28018505, "T", "A", "COSM783", "COSV54042116", 0.094],
#         ["TP53", "S241F", "c.722C>T", "chr17", 7674241, "G", "A", "COSM10812", "COSV52661688", 0.0582],
#     ]
#     known_table = {"data": [], "headers": []}
#     known_table["headers"] = [
#         {"header": "RunID"},
#         {"header": "DNAnr"},
#         {"header": "Gene"},
#         {"header": "Variant AA"},
#         {"header": "CDS mutation"},
#         {"header": "Chr"},
#         {"header": "Pos"},
#         {"header": "Ref"},
#         {"header": "Alt"},
#         {"header": "Legacy ID"},
#         {"header": "Genomic mutation ID"},
#         {"header": "Expected AF"},
#         {"header": "AF"},
#         {"header": "DP"},
#         {"header": "SV Length"},
#     ]
#     for known_line in known:
#         if known_line[1] == "ITD300":
#             for record in pindel_file.fetch(known_line[3], known_line[4] - 1, known_line[4]):
#                 af = int(record.samples.keys()[0].get("AD")[1]) / sum(record.samples.keys()[0].get("AD"))
#                 dp = sum(record.samples.keys()[0].get("AD"))
#                 outline = [sequenceid, sample] + known_line + [float(af), int(dp), record.info["SVLEN"]]
#                 known_table["data"].append(outline)
#                 known_number += 1
#         else:
#             for record in vcf_file.fetch(known_line[3], known_line[4] - 1, known_line[4]):
#                 if record.alts[0] == known_line[6] and record.ref == known_line[5]:
#                     outline = [sequenceid, sample] + known_line + [float(record.info["AF"][0]), int(record.info["DP"])]
#                     known_table["data"].append(outline)
#                     known_number += 1
#
#     return known_table, known_number / len(known)


class VariantTable(object):
    """
    Create Excel-compatible table of variants from VCF.
    Writes SNVs and indels to one table, and SVs (>50 bp) to
    a second one if no path is provided for any VCF with CNV and/or SVs.
    """
    def __init__(self, vcf_snv_indels=None, vcf_cnv_sv=None, min_sv_length=50, min_af=0.005):
        self.vcf_snv_indels = vcf_snv_indels
        self.vcf_cnv_sv = vcf_cnv_sv
        self.min_sv_length = min_sv_length
        self.min_af = min_af
        self.table_snv_indels = {"data": [], "headers": []}
        self.table_cnv_sv = {"data": [], "headers": []}
        self.sv_records = []
        self.sample = ""

    def create_snv_table(self, sequenceid):
        vcf_file = VariantFile(self.vcf_snv_indels)
        self.sample = list(vcf_file.header.samples)[0]
        self.csq_index = index_vep(vcf_file)

        self.table_snv_indels["headers"] = [
            {"header": "Filter"},
            {"header": "RunID"},
            {"header": "DNAnr"},
            {"header": "Gene"},
            {"header": "Chr"},
            {"header": "Pos"},
            {"header": "Ref"},
            {"header": "Alt"},
            {"header": "AF"},
            {"header": "DP"},
            {"header": "Transcript"},
            {"header": "Mutation cds"},
            {"header": "ENSP"},
            {"header": "Consequence"},
            {"header": "COSMIC ids on pos"},
            {"header": "Clinical Significance"},
            {"header": "dbSNP"},
            {"header": "Max Pop AF"},
            {"header": "Max Pop"},
            # {"header": "Artifact Medians (Mutect; Vardict)"},
            # {"header": "Artifact calls (Mutect; Vardict; TotNormals)"},
            # {"header": "Number of SD from background median"},
            # {"header": "Background Median"},
            {"header": "Callers"},
        ]
        for record in vcf_file.fetch():
            record_values = extract_vcf_values(record, self.csq_index)
            var_length = max(len(record.ref), len(record.alts[0]))
            if record_values["af"] > self.min_af and var_length <= self.min_sv_length:
                outline = [
                    record_values["filter_flag"],
                    sequenceid,
                    self.sample,
                    record_values["gene"],
                    record.contig,
                    int(record.pos),
                    record.ref,
                    record.alts[0],
                    record_values["af"],
                    record_values["dp"],
                    record_values["transcript"],
                    record_values["coding_name"],
                    record_values["ensp"],
                    record_values["consequence"],
                    record_values["cosmic"],
                    record_values["clinical"],
                    record_values["rs"],
                    record_values["max_pop_af"],
                    record_values["max_pops"],
                    # record_values["artifact_median"],
                    # record_values["artifact_callers"],
                    # record_values["background_nr_sd"],
                    # record_values["background_median"],
                    record_values["callers"],
                ]
                self.table_snv_indels["data"].append(outline)
            elif record_values["af"] > self.min_af and var_length > self.min_sv_length:
                self.sv_records.append(record)

    def get_table_snv_indels(self):
        return self.table_snv_indels

    def create_sv_table(self, sequenceid):
        if self.vcf_cnv_sv:
            cnv_sv_file = VariantFile(self.vcf_cnv_sv)
            # sample = list(cnv_sv_file.header.samples)[0]
            csq_index = index_vep(cnv_sv_file)
        else:
            cnv_sv_file = None
            sample = None
            csq_index = None

        self.table_cnv_sv["headers"] = [
            {"header": "Filter"},
            {"header": "RunID"},
            {"header": "DNAnr"},
            {"header": "Gene"},
            {"header": "Chr"},
            {"header": "Pos"},
            {"header": "Ref"},
            {"header": "Alt"},
            {"header": "SV length"},
            {"header": "AF"},
            {"header": "DP"},
            {"header": "Transcript"},
            {"header": "Mutation cds"},
            {"header": "ENSP"},
            {"header": "Consequence"},
            {"header": "COSMIC ids on pos"},
            {"header": "Clinical Significance"},
            {"header": "dbSNP"},
            {"header": "Max Pop AF"},
            {"header": "Max Pop"},
            # {"header": "Artifact Medians"},
            # {"header": "Artifact calls (indel; TotNormals)"},
            {"header": "Callers"}
        ]

        if cnv_sv_file:
            for record in cnv_sv_file.fetch():
                self.sv_records.append(record)

        for record in self.sv_records:
            record_values = extract_vcf_values(record, csq_index) if csq_index else extract_vcf_values(record, self.csq_index)
            var_length = max(len(record.ref), len(record.alts[0]))
            if record_values["af"] > self.min_af:
                outline = [
                    record_values["filter_flag"],
                    sequenceid,
                    self.sample,
                    record_values["gene"],
                    record.contig,
                    int(record.pos),
                    record.ref,
                    record.alts[0],
                    var_length,
                    record_values["af"],
                    record_values["dp"],
                    record_values["transcript"],
                    record_values["coding_name"],
                    record_values["ensp"],
                    record_values["consequence"],
                    record_values["cosmic"],
                    record_values["clinical"],
                    record_values["rs"],
                    record_values["max_pop_af"],
                    record_values["max_pops"],
                    # record_values["artifact_median"],
                    # record_values["artifact_callers"],
                    record_values["callers"],
                ]
                self.table_cnv_sv["data"].append(outline)
        else:
            pass

    def get_table_cnv_sv(self):
        return self.table_cnv_sv
