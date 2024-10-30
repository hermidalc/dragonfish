#!/bin/bash

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
while [[ -v CONDA_DEFAULT_ENV ]]; do
    conda deactivate
done
conda activate dragonfish

SNAKEMAKE_CMD="snakemake $@"
echo $SNAKEMAKE_CMD
$SNAKEMAKE_CMD
