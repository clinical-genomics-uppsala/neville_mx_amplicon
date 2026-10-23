# Basecalling

Raw signal (POD5 files) is converted to reads with **[Dorado](https://software-docs.nanoporetech.com/dorado/latest/)** v1.1.1.
The rules used depend on `multisample` (see [Single-sample and multi-sample runs](overview.md#single-sample-and-multi-sample-runs)).
All paths end with the same file:
`basecalling/dorado_duplex/{experiment}_{sample}_{type}_reads.ont_adapt_trim.bam`.

!!! tip
    To skip basecalling and use the reads basecalled by MinKNOW, see [Skipping basecalling](running.md#skipping-basecalling).

## Multiplexed samples (`multisample: True`)

### Duplex basecalling
`basecalling_dorado_duplex_multisamples` runs `dorado duplex` on all the POD5 files of the run, with model
`dna_r10.4.1_e8.2_400bps_sup@v5.2.0` on all GPUs (`--device cuda:all`). Reads are **not** trimmed here, because barcodes are
needed for demultiplexing. The model is copied from `DORADO_MODELS` if available, otherwise downloaded.

### Demultiplexing
`basecalling_dorado_demux` runs `dorado demux` with the sample sheet given in `samplesheet` and kit `SQK-NBD114-24`.
`basecalling_rename_demux_bam` then copies each sample's BAM file to `basecalling/rename_demux_bam/{experiment}_{sample}_{type}_reads.basecalled.bam`.

### Trimming
`basecalling_dorado_trim` runs `dorado trim --sequencing-kit SQK-NBD114-24` to remove adapters and barcodes.

## Single sample (`multisample: False`)

`basecalling_dorado_duplex` runs `dorado duplex` and pipes the reads to `dorado trim`
(`--no-trim-primers --sequencing-kit SQK-LSK114`), which removes adapters but keeps the PCR primers. The primers are soft-clipped
after alignment instead.

`basecalling_dorado_basecaller` (simplex basecalling with `--trim adapters`) is also defined for this mode.

## Conversion to FASTQ

`basecalling_bam_to_fastq` converts the unaligned BAM file to FASTQ with `samtools fastq`, and `basecalling_compress_fastq` gzips it.
This FASTQ file is the input of [Filtlong](prealignment.md).

## Configuration

| Rule | Key | Default | Description |
|------|-----|---------|-------------|
| `basecalling_dorado_duplex_multisamples`, `basecalling_dorado_duplex`, `basecalling_dorado_basecaller` | `model` | `dna_r10.4.1_e8.2_400bps_sup@v5.2.0` | Basecalling model. |
| | `extra` | `--device cuda:all` | Extra Dorado options. |
| | `trim` | see above | Trimming options. |
| `basecalling_dorado_trim` | `trim` | `--sequencing-kit SQK-NBD114-24` | Trimming options after demultiplexing. |

**Cluster resources**: the basecalling rules need a GPU node. At CGU they request `--gres=gres:gpu:4`, 40 threads and 72 hours.

See the [rule list](softwares.md) for the full rule definitions.
