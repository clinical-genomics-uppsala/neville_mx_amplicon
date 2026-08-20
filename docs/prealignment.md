# Pre-alignment Filtering

Before aligning reads to the reference genome, they undergo filtering to discard low-quality data and size-based artifacts.

---

## Tool & Rule

### Filtlong
Filtering is performed by the rule `prealignment_ont_filtlong` using **[Filtlong](https://github.com/rrwick/Filtlong)**. Filtlong filters reads by quality (mean Phred score) and length without requiring a reference genome alignment.

---

## Filtering Rationale

In amplicon sequencing, size filtering is critical to eliminate PCR and sequencing artifacts:

1. **Chimera Prevention**: 
   Very long reads (e.g., > 4000 bp) can represent concatemers/chimeras where multiple amplicons concatenated together during the PCR stages.
2. **Truncated Read Removal**: 
   Very short reads (e.g., < 2000 bp) are often "unfinished" PCR products or reads that got truncated due to early release from the nanopore.
3. **Quality Thresholding**: 
   Filtering out low-quality reads (typically `--min_mean_q 10` or higher) reduces background noise and prevents false-positive variant calling.

---

## Configuration Settings

Parameters are adjustable in `config/config.yaml`:
- `min_length`: Minimal read length (default: 2000 bp)
- `max_length`: Maximal read length (default: 4000 bp)
- `min_mean_q`: Minimum mean read quality score
