"""Annotate VCF variants with artifact statistics from a cohort of chosen samples (ideally normal samples).

Snakemake script. For each variant in a single-sample VCF, look up the matching
entry in an artifact table built from a cohort of chosen samples and prepend
three INFO fields:

- ``Artifact``: number of panel samples in which the variant was observed, one
  value per caller, followed by the panel size (the highest observation count
  in the table).
- ``ArtifactMedian``: median MAF of the variant across the panel, per caller.
- ``ArtifactNrSD``: how many standard deviations the sample AF lies from the
  panel median, per caller. 1000 means the variant was not seen in the panel
  (or the AF is above a median with no spread); 0.0 means the AF is at or below
  a median with no spread.
- If they were only 4 or fewer cohort samples that showed the variant, 
statistics cannot be computed ('no spread' case).

Variants are matched on chromosome, position, ref, alt and variant type. The
specific type (e.g. ``C>T``, ``insertion``, ``deletion``) is tried first, then
the generic type (``SNV`` or ``INDEL``). Variants absent from the table get
zero observations.

Variants at a hotspot position (chromosome + position in the hotspots table)
are never treated as artifacts: they get a ``Hotspot`` INFO entry instead of
the artifact fields. The hotspots table is optional; without it no variant
is treated as a hotspot.

The AF used is read from the FORMAT ``AF`` field of the first sample.

Inputs (snakemake.input):
    vcf: VCF to annotate.
    artifacts: tab-separated artifact table. Header columns are Chromosome,
        pos, variant_type, ref, alt, followed by one (caller, median_MAF,
        sd_MAF) triplet per caller; data rows hold (observations, median, sd).
    hotspots (optional): tab-separated hotspot table whose first two columns are
        chromosome and position.

Output (snakemake.output):
    vcf: annotated VCF with the three artifact INFO fields declared in the
        header.
"""

import pandas as pd
from pysam import VariantFile

def add_artifact_annotation_data(in_vcf_filename, artifacts, hotspots, out_vcf_filename):

    artifact_dict = {}
    hotspots_dict = []
    if hotspots:
        hotspots_dict = [hotspot[0:2] for hotspot in pd.read_table(hotspots).itertuples(index=False, name="Hotspot")]
    # [0:2] picks only "Chromosome", "pos" in the tuple
    header = True
    caller_list = []
    Empty_observation = []
    max_nr_observations = 0
    for line in artifacts:
        lline = line.strip().split("\t")
        if header:
            # Header: Chromosome, pos, variant_type, ref, alt, then one (caller, median_MAF, sd_MAF) triplet per caller
            caller_list = lline[5::3]
            Empty_observation = ["0"] * len(caller_list)
            header = False
            continue
        chrom = lline[0]
        pos = lline[1]
        type = lline[2]
        ref = lline[3]
        alt = lline[4]
        observations = []
        medians = []
        sds = []
        i = 0
        for data in lline[5:]:
            if i % 3 == 0:
                observations.append(data)
            if i % 3 == 1:
                medians.append(data)
            if i % 3 == 2:
                sds.append(data)
            i += 1
        for obs in observations:
            if int(obs) > max_nr_observations:
                max_nr_observations = int(obs)
        artifact_dict[chrom + "_" + pos + "_" + type + "_" + ref + "_" + alt] = [type, ref, alt, observations, medians, sds]

    in_vcf = VariantFile(in_vcf_filename)
    new_header = in_vcf.header
    new_header.info.add(
        "Artifact",
        "1",
        "String",
        "Number of observations of SNV or INDEL in panel samples per caller and finally panel size: %s,Total" %
        ",".join(caller_list)
    )
    new_header.info.add("ArtifactMedian", "1", "String", "Artifact median MAFs in normal panel")
    new_header.info.add("ArtifactNrSD", "1", "String", "Number of Standard Deviations from artifacts in panel median")
    out_vcf = VariantFile(out_vcf_filename, 'w', header=new_header)
    out_vcf.close()
    in_vcf.close()

    out_vcf = open(out_vcf_filename, "a")
    in_vcf = open(in_vcf_filename)
    header = True
    for line in in_vcf:
        if header:
            # out_vcf.write(line)
            if line[:6] == "#CHROM":
                header = False
            continue
        lline = line.strip().split("\t")
        chrom = lline[0]
        pos = lline[1]
        ref = lline[3]
        alt = lline[4]
        if len(ref) == 1 and len(alt) == 1:
            specific_type = f"{ref}>{alt}"
            generic_type = "SNV"
        elif len(alt) > len(ref):
            specific_type = "insertion"
            generic_type = "INDEL"
        else:
            specific_type = "deletion"
            generic_type = "INDEL"
        key = chrom + "_" + pos + "_" + specific_type + "_" + ref + "_" + alt
        if key not in artifact_dict:
            key = chrom + "_" + pos + "_" + generic_type + "_" + ref + "_" + alt
        filter = lline[6]
        format_list = lline[8].split(":")
        format_data = lline[9].split(":")
        AF_index = 0
        i = 0
        for f in format_list:
            if f == "AF":
                AF_index = i
            i += 1
        AF = float(format_data[AF_index])
        Observations = Empty_observation
        Medians = Empty_observation
        SDs = Empty_observation
        if key in artifact_dict:
            Observations = artifact_dict[key][3]
            Medians = artifact_dict[key][4]
            SDs = artifact_dict[key][5]
        INFO = lline[7]

        if (chrom, pos) in hotspots_dict:
            # Never mark as artifacts variants reported in hotspot positions
            Artifact_string = "Hotspot"
        else:
            Artifact_string = "Artifact=" + Observations[0]
            if len(Observations) > 1:
                for obs in Observations[1:]:
                    Artifact_string += "," + obs
                Artifact_string += "," + str(max_nr_observations)
            Artifact_string += ";ArtifactMedian=" + Medians[0]
            if len(Medians) > 1:
                for median in Medians[1:]:
                    Artifact_string += "," + median
            Artifact_string += ";ArtifactNrSD="
            i = 0
            nrsd = 1000
            for sd in SDs:
                if sd == "0":
                    nrsd = 1000
                elif float(sd) == 1000.0 or float(sd) == 0.0:
                    if float(Medians[i]) >= AF:
                        nrsd = 0.0
                    else:
                        nrsd = 1000
                else:
                    nrsd = abs(AF - float(Medians[i])) / float(sd)
                if i == 0:
                    Artifact_string += str(nrsd)
                else:
                    Artifact_string += "," + str(nrsd)
                i += 1
        INFO = Artifact_string + ";" + INFO
        lline[7] = INFO
        out_vcf.write(lline[0])
        for column in lline[1:]:
            out_vcf.write("\t" + column)
        out_vcf.write("\n")
    out_vcf.close()


if __name__ == "__main__":
    log = snakemake.log_fmt_shell(stdout=False, stderr=True)

    add_artifact_annotation_data(
        snakemake.input.vcf,
        open(snakemake.input.artifacts),
        snakemake.input.get("hotspots") or None,
        snakemake.output.vcf,
    )
