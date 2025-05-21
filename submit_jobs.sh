#!/bin/bash
#$ -S /bin/bash
#$ -N db_loader_array
#$ -q d12imppc
#$ -t 1-3
#$ -o logs/job_output.$JOB_ID.$TASK_ID.log
#$ -e logs/job_error.$JOB_ID.$TASK_ID.log
#$ -cwd

PROFILES_DIR="profiles_data"
FILES=($(ls $PROFILES_DIR/*.profiles))
FILE=${FILES[$SGE_TASK_ID - 1]}
GWAS_NAME=$(basename "$FILE" | cut -d"-" -f2)
source /imppc/labs/dnalab/xfarrer/miniforge3/etc/profile.d/conda.sh
conda activate dashboard_project_env

export PYTHONPATH=$PYTHONPATH:/impp/labs/dnalab/share/GEPETO/PolyGenie

mkdir -p logs
python pipeline/process_profile.py "$FILE" "$GWAS_NAME"
