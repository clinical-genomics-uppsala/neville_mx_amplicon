# Variant Annotation & Decomposition

Once variants are called, they undergo normalization and functional annotation to determine their potential clinical impact.

---

## Tools & Rules

### VT (Decomposition)
VCF files containing multi-allelic variants (multiple alternate alleles at the same genomic position) are decomposed using **[VT](https://github.com/atks/vt)**:
- **Rule**: `snv_indels_vt_decompose` (or caller ensembled steps)
- **Method**: Splits multi-allelic records into separate biallelic records. It also left-aligns and normalizes INDELs, ensuring consistent variant representation.

### VEP (Annotation)
Functional annotation is performed using Ensembl's **[Variant Effect Predictor (VEP)](https://www.ensembl.org/info/docs/tools/vep/index.html)**:
- **Rule**: `copy_annotation_vep`
- **Method**: Annotates variant consequences (e.g. missense, stop gained, synonymous), protein changes, HGVS nomenclature, population allele frequencies (e.g., from gnomAD), and clinical significance (ClinVar).
- **Configuration**: Runs in offline mode using a local VEP cache directory (e.g., version 111) to maximize speed and stability on cluster environments.
