# Variant Calling

The pipeline calls single nucleotide variants (SNVs), small insertions/deletions (INDELs), and structural variants (SVs) using four complementary callers.

---

## 1. Variant Callers

### ClairS-TO
- **Purpose**: Call somatic SNVs and INDELs in long-read tumor-only sequencing data.
- **Rule**: `snv_indels_svs_clairs_to`
- **Method**: Evaluates allele frequency distributions against background noise models to identify low-frequency somatic mutations.

### DeepSomatic
- **Purpose**: Deep learning-based somatic variant caller designed for tumor-only/tumor-normal long reads.
- **Rule**: `snv_indels_svs_deepsomatic`
- **Method**: Utilizes convolutional neural networks (CNNs) to call low-frequency somatic variants.

### VarDict
- **Purpose**: Ultra-sensitive variant caller running in amplicon mode.
- **Rule**: `create_background_file_longread_vardict` and ensembled steps.
- **Method**: Sensitive to low-frequency variants, though it may yield higher false-positive rates for short INDELs in homopolymer tracts.

### Sniffles2
- **Purpose**: Structural variant (SV) calling for long-read data.
- **Rule**: `cnv_sv_sniffles2_call`
- **Method**: Detects insertions, deletions, duplications, and inversions.

---

## 2. Variant Allele Frequency (VAF) Calculations

Different callers compute and report VAF (represented as `AF` or `VAF` in output VCFs) using different approaches:

### ClairS-TO
ClairS-TO models somatic and germline variant allele frequencies based on purity and ploidy:
- **Somatic AF**: $\text{AF}_{\text{somatic}} = \frac{p \times V}{p \times C + 2(1-p)}$
- **Germline AF**: $\text{AF}_{\text{germline}} = \frac{p \times V + (1-p)}{p \times C + 2(1-p)}$
  *(where $p$ is tumor purity, $C$ is copy number, and $V$ is variant allele count).*
- In the final VCF format field, `AF` represents the raw calculation: $\text{AF} = \frac{\text{AD}_{\text{alt}}}{\text{DP}}$.

### DeepSomatic
DeepSomatic performs read downsampling during variant candidate evaluation. 
- The reported Allelic Depth (`AD`) and Depth (`DP`) fields in the VCF represent the downsampled reads.
- The reported Allele Frequency is: $\text{AF} = \frac{\text{AD}[1]}{\text{DP}}$.

### VarDict
- VarDict calculates VAF as the ratio of reads supporting the variant to the total coverage at that position.

### Sniffles2
Sniffles2 reports VAF in the INFO field as `VAF`, calculated as:
$$\text{VAF} = \frac{\text{SUPPORT}}{\text{COVERAGE}}$$
- **SUPPORT**: Number of reads directly supporting the structural variant.
- **COVERAGE**: Evaluated dynamically based on the SV type (e.g. coverage at insertion center, averages at start/end for deletions and duplications).
