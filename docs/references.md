# Reference Files

The `neville_mx_amplicon` pipeline requires several reference genomes, gene annotation caches, and systematic noise models to ensure accurate variant calling.

---

## 1. Reference Genome

The primary genomic reference files include:
- **Reference FASTA**: The reference genome sequence (typically GRCh38/hg38). For some steps, a decompressed `.fasta` file is required.
- **Fasta Index (`.fai`)**: Index of the reference genome, created using `samtools faidx`.
- **Sequence Dictionary (`.dict`)**: Sequence dictionary file, created using `samtools dict` or `picard CreateSequenceDictionary`.

---

## 2. VEP Cache

For annotation with the Ensembl Variant Effect Predictor (VEP):
- **VEP Cache Directory**: A local download of the Ensembl VEP database (e.g., version 111) containing transcript coordinates and annotations. This avoids web requests during Snakemake execution.

---

## 3. Systematic Noise Models (Panel of Normals)

Long-read sequencing on Oxford Nanopore platforms generates characteristic sequencing noise (systematic homopolymer errors and mismatch artifacts). The pipeline leverages background files to distinguish real variants from sequencing noise:

### ClairS-TO Background
Somatic variant calling with *ClairS-TO* requires background files generated from normal controls or reference samples. The rule `create_background_file_longread_clairs_to` generates a background frequency model to suppress false positive somatic calls.

### DeepSomatic Background
Somatic calling with *DeepSomatic* utilizes background files generated via `create_background_file_longread_deepsomatic` to model noise parameters in tumor-only settings.

### VarDict Background
*VarDict* runs in amplicon mode and also benefits from background files built from normal reference sequences to filter out low-frequency sequencing errors.
