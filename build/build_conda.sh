#!/usr/bin/env bash
# Language: bash
# Path: `buid/build_conda.sh`
#
# Create a conda environment for the neville_mx_amplicon pipeline.
# Defaults:
#  - ENV_FILE = conda/environment.yml
#  - ENV_NAME = neville_mx_amplicon
#  - RECIPE_DIR = conda/recipe
#  - USE_CONDA_BUILD = false (set to "true" to run conda-build if recipe exists)
#
# Example:
#  bash buid/build_conda.sh --env-file conda/environment.yml --env-name neville_mx_amplicon --use-conda-build true

set -euo pipefail

# Defaults
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
ENV_FILE="${REPO_ROOT}/conda/environment.yml"
ENV_NAME="neville_mx_amplicon"
RECIPE_DIR="${REPO_ROOT}/conda/recipe"
USE_CONDA_BUILD="false"
FORCE="false"
QUIET="false"

print_usage() {
  cat <<EOF
Usage: bash buid/build_conda.sh [options]

Options:
  --env-file PATH       Path to environment yaml (default: ${ENV_FILE})
  --env-name NAME       Name of conda environment to create (default: ${ENV_NAME})
  --recipe-dir PATH     Path to conda recipe directory (default: ${RECIPE_DIR})
  --use-conda-build VAL If "true" and recipe exists, run conda-build (default: ${USE_CONDA_BUILD})
  --force               Overwrite existing environment without prompt
  --quiet               Minimise output
  -h|--help             Show this help
EOF
}

# Simple arg parser
while [[ $# -gt 0 ]]; do
  case "$1" in
    --env-file) ENV_FILE="$2"; shift 2;;
    --env-name) ENV_NAME="$2"; shift 2;;
    --recipe-dir) RECIPE_DIR="$2"; shift 2;;
    --use-conda-build) USE_CONDA_BUILD="$2"; shift 2;;
    --force) FORCE="true"; shift;;
    --quiet) QUIET="true"; shift;;
    -h|--help) print_usage; exit 0;;
    *) echo "Unknown option: $1"; print_usage; exit 2;;
  esac
done

info() {
  if [ "${QUIET}" != "true" ]; then
    echo "[INFO] $*"
  fi
}

warn() { echo "[WARN] $*"; }
error() { echo "[ERROR] $*" >&2; exit 1; }

# Check env file exists
if [ ! -f "${ENV_FILE}" ]; then
  error "Environment file not found: ${ENV_FILE}"
fi

# Detect mamba/conda
CONDA_CMD=""
if command -v mamba >/dev/null 2>&1; then
  CONDA_CMD="mamba"
elif command -v conda >/dev/null 2>&1; then
  CONDA_CMD="conda"
else
  error "Neither mamba nor conda found in PATH. Install one to proceed."
fi

info "Using ${CONDA_CMD} to create environment"
info "Environment file: ${ENV_FILE}"
info "Environment name: ${ENV_NAME}"

# If environment already exists
if "${CONDA_CMD}" env list | awk '{print $1}' | grep -qx "${ENV_NAME}"; then
  if [ "${FORCE}" = "true" ]; then
    info "Removing existing environment ${ENV_NAME} (forced)"
    "${CONDA_CMD}" env remove -n "${ENV_NAME}" -y
  else
    read -r -p "Environment ${ENV_NAME} exists. Remove and recreate? [y/N]: " ans
    case "${ans}" in
      [Yy]*) "${CONDA_CMD}" env remove -n "${ENV_NAME}" -y;;
      *) info "Aborting."; exit 0;;
    esac
  fi
fi

# Create the environment
info "Creating environment ${ENV_NAME} from ${ENV_FILE}"
if [ "${CONDA_CMD}" = "mamba" ]; then
  mamba env create -f "${ENV_FILE}" -n "${ENV_NAME}" --yes
else
  conda env create -f "${ENV_FILE}" -n "${ENV_NAME}" --yes
fi

# Export a pinned environment file / lock
EXPORT_FILE="${REPO_ROOT}/conda/env_export.yml"
info "Exporting resolved environment to ${EXPORT_FILE}"
# Activate then export to ensure exact package builds if conda is used
# Use conda run to avoid activating in current shell
if command -v conda >/dev/null 2>&1; then
  conda run -n "${ENV_NAME}" --no-capture-output python - <<PYTHON 2>/dev/null || true
# dummy to ensure environment is writable by conda run
PYTHON
  conda env export -n "${ENV_NAME}" --no-builds > "${EXPORT_FILE}"
else
  # fallback: use mamba/conda env export via runscript
  "${CONDA_CMD}" env export -n "${ENV_NAME}" --no-builds > "${EXPORT_FILE}"
fi

info "Export saved. You can use ${EXPORT_FILE} to reproduce exact package versions."

# Optionally build conda package if recipe present and requested
if [ "${USE_CONDA_BUILD}" = "true" ] && [ -d "${RECIPE_DIR}" ]; then
  if ! command -v conda-build >/dev/null 2>&1; then
    warn "conda-build not found. Skipping conda-build step. Install conda-build to enable packaging."
  else
    info "Running conda-build on recipe at ${RECIPE_DIR}"
    # prefer conda-forge and bioconda channels commonly used in bioinformatics pipelines
    conda-build "${RECIPE_DIR}" -c conda-forge -c bioconda -c defaults
    info "conda-build finished. Built packages are stored in the conda-bld directory."
  fi
else
  if [ "${USE_CONDA_BUILD}" = "true" ]; then
    warn "Requested conda-build but recipe directory not found: ${RECIPE_DIR}. Skipping."
  fi
fi

info "Done. Activate the environment with:"
echo
echo "  conda activate ${ENV_NAME}"
echo

exit 0