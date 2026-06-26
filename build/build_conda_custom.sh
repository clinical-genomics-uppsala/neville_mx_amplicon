#!/usr/bin/env bash
set -e

eval "$(conda shell.bash hook)"

PIPELINE_NAME="neville_mx_amplicon"
PIPELINE_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/neville_mx_amplicon.git"
TAG_OR_BRANCH="offline-pack"
# CONFIG_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/neville_mx_amplicon.git"
# CONFIG_VERSION="develop"
PYTHON_VERSION="3.9"

# Clone git of neville_mx_amplicon to configure conda environment
git clone --branch ${TAG_OR_BRANCH} ${PIPELINE_GITHUB_REPO}
cd ${PIPELINE_NAME}

# Create and activate conda environment in the current directory, then install pipeline requirements
conda create --prefix ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env python=${PYTHON_VERSION} -y
conda activate ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env
conda install -c conda-forge pip -y

if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH} ];
then
    rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}
fi

# The directory ${PIPELINE_NAME}_${TAG_OR_BRANCH} is created for the files that are to be packaged and transferred
# elsewhere:
# - the pipeline code for neville_mx_amplicon
# - the conda environment used to run the pipeline
# - the snakemake-wrappers
# - the hydra-genetics modules that are used in the pipeline: snv_indels, cnv_sv, annotation, filtering, qc, references.

mkdir -p ${PIPELINE_NAME}_${TAG_OR_BRANCH}

# Clone git of neville_mx_amplicon
git clone --branch ${TAG_OR_BRANCH} ${PIPELINE_GITHUB_REPO} ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}
./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env/bin/pip3 install -r ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/requirements.txt
conda pack --prefix ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env -o ${PIPELINE_NAME}_${TAG_OR_BRANCH}/env.tar.gz

# Clone snakemake-wrappers and hydra-genetics
mkdir -p ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics

# Clone wrappers
git clone https://github.com/snakemake/snakemake-wrappers.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/snakemake-wrappers

# Clone hydra modules
git clone https://github.com/hydra-genetics/annotation.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/annotation
git clone https://github.com/hydra-genetics/cnv_sv.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/cnv_sv
git clone https://github.com/hydra-genetics/qc.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/qc
git clone https://github.com/hydra-genetics/snv_indels.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/snv_indels
git clone https://github.com/hydra-genetics/references.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/references

## Save DockerHub paths in YAML to later create Singularity images
mv ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/config.yaml \
${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/containers.yaml

## Download the config files from the config repo
# git clone --branch ${CONFIG_VERSION} ${CONFIG_GITHUB_REPO} neville_mx_amplicon_config
# TODO

# Pack all cloned repositories
tar -zcvf ${PIPELINE_NAME}_${TAG_OR_BRANCH}.tar.gz ${PIPELINE_NAME}_${TAG_OR_BRANCH}

# Download containers
conda activate ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env
hydra-genetics prepare-environment create-singularity-files -c ./${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/containers.yaml -o apptainer_cache

# Download references
for reference_config in "$@"
do
    hydra-genetics --debug references download -o design_and_ref_files -v $reference_config
done

conda deactivate

if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH}_env ];
then
    rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}_env
fi

if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH} ];
then
    rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}
fi
