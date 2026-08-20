# Frequently Asked Questions

---

## 1. How do I bypass Dorado basecalling if I already have basecalled BAM files?
If basecalling was already completed on the sequencing computer by MinKNOW:
1. Merge the individual pass BAMs in the `bam_pass/` directory into a single BAM file using `samtools merge`.
2. Place the merged BAM file inside the path designated in the `units.tsv` file under the `bam` column.
3. Adjust the basecalling rule inputs in `workflow/rules/basecalling.smk` to point directly to these files.

---

## 2. Why does Dorado fail with "Failed to load NVML / no CUDA devices available"?
This occurs when Dorado is requested to run on a GPU but can't find or access CUDA:
- Ensure you have requested GPU resources in the Slurm profile (e.g. `--gres=gpu:4` or similar).
- Verify that Singularity/Apptainer is run with the `--nv` option to expose host NVIDIA drivers to the container.

---

## 3. How does the pipeline handle coverage for overlapping TP53 amplicons?
Amplicons like `TP53_D` and `TP53_J` overlap significantly. Standard coverage calculations would double-count reads in the overlapping segments. 
- The pipeline utilizes custom BED files defining unique non-overlapping segments.
- For complete overlaps (e.g. `D` and `J`), the coverage for `J` is calculated in post-processing by subtracting the coverage of `D` from the total coverage of the combined overlapping region.
