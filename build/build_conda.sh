#!/usr/bin/env bash
set -e

eval "$(conda shell.bash hook)"

PIPELINE_NAME="neville_mx_amplicon"
PIPELINE_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/neville_mx_amplicon.git"
TAG_OR_BRANCH="offline-pack"
# CONFIG_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/neville_mx_amplicon.git"
# CONFIG_VERSION="develop"
PYTHON_VERSION="3.9"
PACK_NAME="${PIPELINE_NAME}_packs"

# Choose what build steps to execute
BUILD_CONDA_ENV=true
BUILD_APPTAINERS=false
BUILD_SMK_WRAPPERS=true
BUILD_HYDRA_MODULES=true
BUILD_REFERENCES=true

mkdir -p ${PACK_NAME}
cd ${PACK_NAME}

# Create and activate conda environment in the current directory, then install pipeline requirements
if ${BUILD_CONDA_ENV};
then
	conda create --prefix ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env python=${PYTHON_VERSION} -y
	conda activate ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env
	conda install -c conda-forge pip -y
fi

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
# tmp directory required to run some tools in the pipeline
mkdir -p ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/tmp


if ${BUILD_CONDA_ENV};
then
	./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env/bin/pip3 install -r ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/requirements.txt
	conda pack --prefix ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env -o ${PIPELINE_NAME}_${TAG_OR_BRANCH}/env.tar.gz
fi

# Clone snakemake-wrappers and hydra-genetics
mkdir -p ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics

# Clone wrappers
if ${BUILD_SMK_WRAPPERS};
then
	git clone https://github.com/snakemake/snakemake-wrappers.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/snakemake-wrappers
fi

# Clone hydra modules
if ${BUILD_HYDRA_MODULES};
then
	git clone https://github.com/hydra-genetics/annotation.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/annotation
	git clone https://github.com/hydra-genetics/cnv_sv.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/cnv_sv
	git clone https://github.com/hydra-genetics/filtering.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/filtering
	git clone https://github.com/hydra-genetics/qc.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/qc
	git clone https://github.com/hydra-genetics/snv_indels.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/snv_indels
	git clone https://github.com/hydra-genetics/references.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/references
fi

## Save DockerHub paths in YAML to later create Singularity images
mv ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/config.yaml \
${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/containers.yaml

## Download the config files from the config repo
# git clone --branch ${CONFIG_VERSION} ${CONFIG_GITHUB_REPO} neville_mx_amplicon_config
# TODO
sed -i -E "s/TAG_OR_BRANCH/${TAG_OR_BRANCH}/g" ./${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/profiles/miarka/config.yaml
sed -i -E "s/TAG_OR_BRANCH/${TAG_OR_BRANCH}/g" ./${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/profiles/slurm/config.yaml

# Pack all cloned repositories
tar -zcvf ${PIPELINE_NAME}_${TAG_OR_BRANCH}.tar.gz ${PIPELINE_NAME}_${TAG_OR_BRANCH}

# Download containers
if ${BUILD_APPTAINERS};
then
	conda activate ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env
	hydra-genetics prepare-environment create-singularity-files -c ./${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/config/containers.yaml -o apptainer_cache
fi

# Download references
if ${BUILD_REFERENCES};
then
	for reference_config in "$@"
	do
		hydra-genetics --debug references download -o design_and_ref_files -v $reference_config
	done
fi

conda deactivate

if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH}_env ];
then
    rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}_env
fi

if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH} ];
then
    rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}
fi

# Compress data
tar -czvf design_and_ref_files.tar.gz design_and_ref_files