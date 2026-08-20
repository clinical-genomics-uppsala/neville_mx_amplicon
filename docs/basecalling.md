# Basecalling & Demultiplexing

The pipeline starts by basecalling raw Nanopore electrical signal data and separating barcoded runs into individual sample streams.

---

## Tools & Rules

### Dorado Basecaller
Basecalling is executed using ONT's official basecaller **[Dorado](https://github.com/nanoporetech/dorado)**. It processes raw `.pod5` files and generates unmapped basecalled reads in BAM format.
- **Simplex basecalling**: Handled by the rule `basecalling_dorado_basecaller`.
- **Duplex basecalling**: Handled by the rules `basecalling_dorado_duplex` and `basecalling_dorado_duplex_multisamples` to call high-accuracy duplex reads (sequenced from both strands).

### Dorado Demux
For barcoded runs, the rule `basecalling_dorado_demux` separates the merged basecalled BAM file into barcoded files based on a user-provided sample sheet.

### Dorado Trim
Adapters and barcodes are clipped from reads using `basecalling_dorado_trim` to ensure they do not interfere with genomic alignment.

### Bam to Fastq Conversion
The rule `basecalling_bam_to_fastq` converts the basecalled/trimmed BAM files into compressed FASTQ files using `samtools fasta` so they can be processed by Filtlong.

---

## Skipping Basecalling
If MinKNOW has already performed basecalling and quality-filtering on the sequencer, you can bypass this step. 
- Merge the BAM files from the run's `bam_pass` directory using `samtools merge`.
- Edit the basecalling rules in `workflow/rules/basecalling.smk` to point the inputs directly to the pre-existing BAM file.
