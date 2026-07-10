#!/bin/bash
#SBATCH -A ngi2024001
#SBATCH -p gpu
#SBATCH --gres=gpu:a100:2
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=00:05:00
#SBATCH -J gpu_diag
#SBATCH -o gpu_diag.%j.out
#SBATCH -e gpu_diag.%j.err

echo "===== HOST ====="
hostname
date

echo
echo "===== SLURM ====="
env | grep '^SLURM' | sort

echo
echo "===== CUDA ====="
echo "CUDA_VISIBLE_DEVICES=${CUDA_VISIBLE_DEVICES:-<unset>}"
echo "NVIDIA_VISIBLE_DEVICES=${NVIDIA_VISIBLE_DEVICES:-<unset>}"
echo "SLURM_JOB_GPUS=${SLURM_JOB_GPUS:-<unset>}}
echo "SLURM_STEP_GPUS=${SLURM_STEP_GPUS:-<unset>}
echo "SLURM_GPUS=${SLURM_GPUS:-<unset>}}

echo
echo "===== NVIDIA-SMI ====="
nvidia-smi

echo
echo "===== APPTAINER ====="

singularity exec --nv \
/proj/ngi2024001/nobackup/bin/wp2_mxa/apptainer_cache/ontresearch__dorado_latest.sif \
bash -c '

echo "Container hostname:"
hostname
env | sort > gpu_diagnostics.env

echo
echo "Environment:"
env | grep -E "CUDA|NVIDIA|SLURM" | sort

echo
echo "GPU devices:"
ls -l /dev/nvidia*

echo
echo "Driver:"
cat /proc/driver/nvidia/version

echo
echo "Dorado:"
dorado --version

echo
echo "Testing CUDA initialization..."

dorado basecaller \
    --device cuda:all \
    /proj/ngi2024001/nobackup/projects/wp3/GMS_str/dorado_models/dna_r10.4.1_e8.2_400bps_sup@v5.2.0 \
    /proj/ngi2024001/nobackup/inbox/wp2_mxa/CGU_2024_05_M32_260421_JNZ/M32/20260421_1315_MN48987_FBE71798_59691875/pod5 \
    --emit-fastq \
    >/dev/null
'