# Alignment & Primer Soft-clipping

Filtered reads are aligned to the reference genome, sorted, indexed, and processed to remove primer binding sequences.

---

## Tools & Rules

### Dorado Align
Alignment is performed by the rule `alignment_ont_dorado_align` using **[Dorado](https://github.com/nanoporetech/dorado)** (which wraps `minimap2` internally for long-read alignment). Read groups (@RG tags) are populated dynamically using run metadata from `units.tsv`.

### Samtools Sort & Index
Aligned BAM files are coordinate-sorted with `alignment_ont_bam_sort` and indexed with `alignment_ont_bam_index` using **[Samtools](http://www.htslib.org/)**.

### Samtools ampliconclip (Soft-clipping)
The rule `alignment_ont_bam_softclip` uses `samtools ampliconclip` to soft-clip primer coordinates from the aligned reads. 

- **Why clip primers?** 
  Since primers are synthetic oligonucleotides that bind to targets, they do not represent the native patient sequence. If mutations occur within the primer binding region of the template DNA, the sequencing reads will show the wild-type primer sequence instead of the mutation, leading to false negatives. Soft-clipping masks these primer regions from variant callers.
- **Reference**: Requires a primer BED file specifying the forward and reverse primer coordinates for each amplicon.

---

## Downstream Processing

After clipping, BAM files are indexed, and the rule `alignment_ont_split_bam_by_target` splits the alignments into target-specific regions. This enables faster, localized execution for variant calling and coverage verification.
