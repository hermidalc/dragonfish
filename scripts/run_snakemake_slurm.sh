#!/bin/bash

CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
while [[ -v CONDA_DEFAULT_ENV ]]; do
    conda deactivate
done

i=1
SNAKEMAKE_OPTS=()
while [[ $i -le $# ]]; do
    if [[ ${!i} == "--sbatch-opts="* ]]; then
        SBATCH_OPTS=${!i#*=}
    elif [[ ${!i} == "--sbatch-opts" ]]; then
        i=$((i + 1))
        SBATCH_OPTS=${!i}
    else
        SNAKEMAKE_OPTS+=(${!i})
    fi
    i=$((i + 1))
done

SCRIPT_DIR=$(dirname $(realpath -s $0))
PROJECT_DIR=$(realpath $SCRIPT_DIR/../)

mkdir -p $PROJECT_DIR/logs

SBATCH_CMD="sbatch \
--chdir=$PROJECT_DIR \
--output=$PROJECT_DIR/logs/run_snakemake_%j.log \
$SBATCH_OPTS \
$SCRIPT_DIR/run_snakemake.sh ${SNAKEMAKE_OPTS[@]}"
echo $SBATCH_CMD
$SBATCH_CMD
