# Quality Control & Reports

The pipeline includes a comprehensive quality control suite to monitor library yield, amplicon coverage, and sequencing quality.

---

## 1. Quality Control Metrics & Tools

### Read Quality (Sequali)
- **Tool**: **[Sequali](https://github.com/rhpvorderman/sequali)**
- **Rule**: `qc_ont_sequali`
- **Description**: Analyzes read length, quality score distributions, and nucleotide compositions, specifically optimized for ONT long-read sequencing.

### Coverage Assessment (Mosdepth)
- **Tool**: **[Mosdepth](https://github.com/brentp/mosdepth)**
- **Rule**: `qc_ont_mosdepth_overlap`
- **Description**: Computes sequencing depth across specified target regions. Because amplicons can overlap (especially in genes like *TP53*), custom BED files covering *only* unique, non-overlapping amplicon segments are used to estimate the coverage of each specific target accurately.

### Off-target Rate (Picard)
- **Tool**: **[Picard](https://broadinstitute.github.io/picard/)**
- **Rule**: `qc_ont_picard_bed_to_interval_list` and downstream hs metrics rules.
- **Description**: Evaluates the percentage of reads mapping off-target, showing how efficiently the multiplexed PCR amplification targeted the desired loci.

---

## 2. Yield and Pool Balance Assessment

In a multiplexed PCR setup, checking the balance of reads is crucial:
- **Within-pool Balance**: For a pool containing 5 amplicons, each amplicon should ideally capture about 20% of the reads. If one amplicon dominates, it suggests PCR bias.
- **Between-pool Balance**: Across 3 pools, each pool should ideally represent 1/3 (~33%) of the total sequencing yield.

The pipeline computes these distributions using custom scripts and writes them out as CSV/PNG files (e.g. `timestep_coverage_per_amplicon` plots) to help fine-tune pipetting quantities for future library preparations.

---

## 3. Compiled Reports

All individual reports are compiled into two final deliverables:
1. **MultiQC Report**: An interactive HTML page combining stats from Sequali, Picard, and Mosdepth.
2. **Excel Spreadsheet (`{sample}_{type}.xlsx`)**: A clinical-friendly spreadsheet containing prefiltered variant calls, per-amplicon coverage values, and pool yield statistics.
