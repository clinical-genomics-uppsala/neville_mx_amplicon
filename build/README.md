Build for packaging, the script must be run from a directory where the pipeline was fully cloned.

```bash
bash build/build_conda.sh ../config/references/references.hg38.md5sums.yaml
```

Expected file structure while building:

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
│   │   └── neville_mx_amplicon
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

Expected file structure after building completed:

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