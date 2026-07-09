# Multi-stage build for packaging

## Preparation

For improved flexibility, you may choose what parts of the pipeline must be (re)packaged.
Set the following bash variables in `build/build_conda.sh` to true or false depending on the stages that you want to build:

```
BUILD_CONDA_ENV=<true|false>
BUILD_APPTAINERS=<true|false>
BUILD_SMK_WRAPPERS=<true|false>
BUILD_HYDRA_MODULES=<true|false>
BUILD_REFERENCES=<true|false>
```

For instance, this will (re)build all stages:
```
BUILD_CONDA_ENV=true
BUILD_APPTAINERS=true
BUILD_SMK_WRAPPERS=true
BUILD_HYDRA_MODULES=true
BUILD_REFERENCES=true
```

But this will rebuild only the code for the pipeline with the suitable smk_wrappers and hydra-genetics modules:

```
BUILD_CONDA_ENV=false
BUILD_APPTAINERS=false
BUILD_SMK_WRAPPERS=true
BUILD_HYDRA_MODULES=true
BUILD_REFERENCES=false
```

## Run the build script for the selected stages

The build script must be run from a directory where the pipeline was fully cloned:

```bash
bash build/build_conda.sh ../config/references/references.hg38.md5sums.yaml
```

Expected file structure while building all stages:

```
neville_mx_amplicon/
├── ...
│   ├── ...
├── neville_mx_amplicon_packs/ 
│   ├── apptainer_cache/
│   ├── design_and_ref_files/
│   │   └── neville_mx_amplicon
│   ├── neville_mx_amplicon_offline-pack/
│   │   ├── hydra-genetics
│	│	├── neville_mx_amplicon
│   │   └── snakemake-wrappers
│   ├── neville_mx_amplicon_offline-pack_env/
│   │   ├── bin
│   │   ├── compiler_compat
│   │   ├── conda-meta
│   │   ├── etc
│   │   ├── include
│   │   ├── lib
│   │   ├── share
│   │   ├── ssl
│   │   └── x86_64-conda-linux-gnu  
│   └── neville_mx_amplicon_offline-pack.tar.gz 
├── ...
```

Expected file structure after building all stages completed:

```
neville_mx_amplicon/
├── ...
│   ├── ...
├── neville_mx_amplicon_packs/ 
│   ├── apptainer_cache/
│   ├── design_and_ref_files/
│   │   └── neville_mx_amplicon  
│   └── neville_mx_amplicon_offline-pack.tar.gz 
├── ...
```
