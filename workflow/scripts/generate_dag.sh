# Usage: `cd neville_mx_amplicon && bash workflow/scripts/generate_dag.sh`
snakemake -s workflow/Snakefile --configfile config/config.yaml --dag | dot -Tpdf > docs/images/dag.pdf
snakemake -s workflow/Snakefile --configfile config/config.yaml --dag | dot -Tpng > docs/images/dag.png
snakemake -s workflow/Snakefile --configfile config/config.yaml --dag > docs/images/dag.dot

snakemake -s workflow/Snakefile --configfile config/config.yaml --rulegraph | dot -Tpdf > docs/images/rulegraph.pdf
snakemake -s workflow/Snakefile --configfile config/config.yaml --rulegraph | dot -Tpng > docs/images/rulegraph.png
snakemake -s workflow/Snakefile --configfile config/config.yaml --rulegraph > docs/images/rulegraph.dot