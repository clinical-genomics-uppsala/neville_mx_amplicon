# Running the pipeline

## Before you start

1. Install the pipeline and edit `config/config.yaml` (see [Pipeline setup and configuration](setup.md)).
2. Create `samples.tsv` and `units.tsv` in the working directory (see [Sample files](setup.md#sample-files)).
3. Create a `tmp/` folder in the working directory.
4. Choose or adapt a [profile](setup.md#profiles).

## Expected run folder layout

The pipeline looks for the POD5 files in `{runfolder}/{batchid}/{runid}/pod5`, which matches the MinKNOW output:

```
<runfolder>/
└── <batchid>/
    └── <runid>/                      # e.g. 20251017_1440_MN48987_FBD39643_37b961fb
        ├── pod5/
        ├── bam_pass/<barcode>/*.bam
        ├── sequencing_summary_*.txt  # used by pycoQC
        └── Samplesheet_<batchid>.csv # multiplexed runs only
```

## Command

```bash
source .venv/bin/activate
snakemake --profile profiles/<profile>/ -s workflow/Snakefile \
    --configfile config/config.yaml \
    --config runfolder=<runfolder> batchid=<batchid> runid=<runid> multisample=True samplesheet=<samplesheet.csv>
```

Useful options:

| Option | Effect |
|--------|--------|
| `-n` | Dry run: list the jobs without running them. |
| `--notemp` | Keep all intermediate files. |
| `--use-singularity --singularity-args "..."` | Needed when not set in the profile. |
| `--configfiles a.yaml b.yaml` | Merge several config files; later files override earlier ones. |

## Multiplexed barcoded samples (MinION)

`workflow/scripts/start_marvin_multisample.sh` wraps all the steps for a MinION run with barcoded samples on the CGU cluster:

```bash
bash workflow/scripts/start_marvin_multisample.sh <runFolder> <batchId> <flowcellId>
```

The script:

1. finds the run folder whose name contains `<flowcellId>` under `<runFolder>/<batchId>/`;
2. reads the sample sheet `Samplesheet_<batchId>.csv` from that folder;
3. for each sample, merges the MinKNOW BAM files of its barcode and runs `hydra-genetics create-input-files` on the merged file;
4. combines the per-sample files into one `samples.tsv` and one `units.tsv`;
5. starts Snakemake with `multisample=True`, which basecalls all POD5 files with Dorado and demultiplexes them with the sample sheet.

!!! note
    Dorado needs a [user-defined sample sheet](https://software-docs.nanoporetech.com/dorado/latest/barcoding/sample_sheet/)
    to demultiplex. The script expects the columns
    `position_id,flow_cell_id,kit,experiment_id,sample_id,alias,barcode`, with a header line.

`start_marvin_multisample_from_fastq.sh` and `start_wp2_mxa_minion_miarka.sh` are variants of the same script for other set-ups.

## Single sample (Flongle)

`workflow/scripts/start_marvin_single.sh` runs one non-barcoded sample with `multisample=False`.
ONT is discontinuing Flongle flowcells, so this script is kept for legacy runs only.

## Skipping basecalling

If the MinKNOW basecalling settings are good enough, you can start from the BAM files in `bam_pass/`:

1. merge them with `samtools merge -o merged.bam bam_pass/*.bam`;
2. copy the result to `basecalling/dorado_duplex/{batchid}_{sample}_T_reads.ont_adapt_trim.bam` in the working directory.

Snakemake then starts from this file and does not run the basecalling rules.

## GPU basecalling

Dorado runs on GPU. When the pipeline basecalls on a cluster:

* request GPUs for the basecalling rules (`gres` and GPU `partition` in `resources.yaml`, or `set-resources` in a Slurm profile);
* pass `--nv` to Singularity/Apptainer.

Without `--nv` Dorado fails with:

```
[error] device string set to cuda:all but no CUDA devices available.
```

The message `INFO: Could not find any nv files on this host!` on CPU-only rules is harmless.
`workflow/scripts/gpu_diagnostics.sh` helps check that the GPU is visible from the container.
