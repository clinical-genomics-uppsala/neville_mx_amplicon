import snakemake
from heapq import merge
from pybedtools import BedTool


def non_overlapping(data):
    out = []
    starts = sorted([(i[0], 1) for i in data])  # start of interval adds a layer of overlap
    ends = sorted([(i[1], -1) for i in data])   # end removes one
    layers = 0
    current = []
    for value, event in merge(starts, ends):    # sorted by value, then ends (-1) before starts (1)
        layers += event
        if layers == 1:  # start of a new non-overlapping interval
            current.append(value + 1)
        elif current:  # we either got out of an interval, or started an overlap
            current.append(value - 1)
            out.append(current)
            current = []
    return out

def bed_to_list(bed_path: str, chrom) -> list:
    intervals = []
    with open(bed_path, "r") as bed_file:
        for line in bed_file:
            if line.startswith(f"chr{chrom}"):
                intervals.append(list(map(int, line.rstrip("\n").split("\t")[1:3])))
    return intervals


def list_to_bed(intervals, chrom, bedout, first='A'):
    with open(bedout, "w") as bed_file:
        for i, interval in enumerate(intervals):
            bed_file.write("\t".join([f"chr{chrom}"] + list(map(str, interval)) + [chr(ord(first) + i)]) + "\n")


bedHI = [[7686781, 7689769], [7688621, 7691634]]


if __name__ == "__main__":
    #print_reduced_intervals("/home/camille/ampliconthemato/amplic_ont_hemato/data/primer_data/amplicons.bed")
    spans = bed_to_list("/home/camille/ampliconthemato/amplic_ont_hemato/data/primer_data/amplicons.bed", 17)
    print(spans)
    print(non_overlapping(spans))
    list_to_bed(non_overlapping(spans), 17,"/home/camille/ampliconthemato/amplic_ont_hemato/data/primer_data/TP53_nonoverlapping.bed")

