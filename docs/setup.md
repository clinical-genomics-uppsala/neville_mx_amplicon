# Pipeline Setup and Configuration

Before executing the pipeline, the environment must be configured, sample files prepared, and design coordinates defined in the configuration.

---

## 1. Virtual Environment Setup

Create a virtual environment and install the required dependencies:
```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -I -r requirements.txt
```
Additionally, create a `tmp` directory at the root of the repository to store temporary files:
```bash
mkdir -p tmp
```

---

## 2. Sample Data Files (`samples.tsv` and `units.tsv`)

The pipeline relies on `samples.tsv` and `units.tsv` configuration files to resolve inputs. 
You can use the helper command `hydra-genetics create-input-files` to initialize these files.

### samples.tsv
A single-column table listing the unique sample/patient IDs.

| sample |
|--------|
| D25-test007 |

### units.tsv
Provides sequencing metadata and pointers to the raw unmapped BAM files.

| Column | Description | Example |
|--------|-------------|---------|
| `sample` | The unique sample/patient ID matching `samples.tsv` | `D25-test007` |
| `type` | Data type identifier (T for tumor) | `T` |
| `platform` | Sequencing platform (ONT) | `ONT` |
| `machine` | Specific sequencer model or identifier | `MinION` |
| `processing_unit` | Flowcell identifier | `FAW12345` |
| `run_id` | Unique sequencing run ID | `run_abc123` |
| `barcode` | Barcode identifier or `NA` if single sample | `NA` |
| `methylation` | Use methylation information (`Yes`/`No`) | `No` |
| `basecalling_model` | Basecalling model name used in MinKNOW | `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` |
| `bam` | Absolute path to the unmapped/basecalled BAM file | `/path/to/reads.bam` |

---

## 3. BED Design Files

The multiplexed amplicon setup requires several design-specific BED files to instruct rules on clipping and coverage regions.

1. **Amplicon BED File**: 
   Defines the genomic region of each amplicon. Coordinates must be sorted chromosomally, with no empty lines.
2. **Primer BED File**: 
   Contains the coordinates of the forward and reverse primers for each amplicon. Required for soft-clipping with *Samtools*.
3. **Overlapping Region BED Files**: 
   For TP53 (or other overlapping designs), specific coordinates covering *only* the unique non-overlapping segments of amplicons are needed to ensure mosdepth accurately calculates per-amplicon coverage.
4. **Caller-specific BED Files**: 
   Restricts variant calling to desired areas to increase performance for *ClairS-TO*, *DeepSomatic*, and *VarDict*.

Examples of these files are located under `.tests/integration/test_data/bedfiles/`.

---

## 4. Configuration yaml (`config/config.yaml`)

The main configuration file coordinates inputs, reference genomes, parameters, and outputs. A template config is provided in `config/config.yaml`.

Key sections include:
- `reference`: Paths to reference fasta, indexes, VEP caches.
- `bedfiles`: Paths to design BED files.
- `parameters`: Tuning options for Filtlong (read length, quality), variant callers, and annotation.

---

## 5. Profile Configurations

Snakemake profiles are stored in the `profiles/` directory:
- `profiles/local/`: For running locally on a single machine or laptop.
- `profiles/slurm/`: For running on HPC clusters using SLURM scheduler resources.
