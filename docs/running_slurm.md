# Running on HPC (Slurm)

In clinical diagnostics and production environments, the pipeline is typically run on an HPC cluster managed by a Slurm workload scheduler.

---

## 1. Slurm Profile Configuration

Modify `profiles/slurm/config.yaml` to specify default cluster options, partition parameters, and thread limits suitable for your environment.

### GPU Allocation for Basecalling
If you are performing basecalling with *Dorado* on GPU nodes, you must define the correct resources to request the GPUs (e.g., passing `--gres=gres:gpu:4` to DRMAA). 

Additionally, pass the `--nv` option to Singularity/Apptainer to mount the NVIDIA GPU drivers inside the container:
```yaml
singularity-args: "--nv -B /tmp -B /scratch"
```
*If `--nv` is omitted, Dorado will fail with CUDA device connection errors.*

---

## 2. Cluster Submission Scripts

Two example bash scripts are provided in `workflow/scripts/` to automate the configuration, environment loading, input merging, and Snakemake execution:

### Multisample Runs (`start_marvin_multisample.sh`)
Used for sequencing runs on MinION flowcells containing multiplexed barcoded samples.
- **Requirements**: Requires a [Dorado sample sheet](https://software-docs.nanoporetech.com/dorado/latest/barcoding/sample_sheet/) defining barcodes to sample names.
- **Workflow**: Performs basecalling, demultiplexing, and outputs barcoded BAMs before running the rest of the pipeline.

### Single Sample Runs (`start_marvin_single.sh`)
Used for Flongle flowcells or single non-barcoded samples.
- **Note**: Flongle flowcells are deprecated by ONT. This script is maintained primarily for legacy and archival purposes.
