#!/bin/bash

set -e  # Exit immediately if a command exits with a non-zero status

polygenie_dir='/imppc/labs/dnalab/share/GEPETO/PolyGenie'
genopred_dir='/imppc/labs/dnalab/share/software/GenoPred'

# Get the number of cores from the first argument
NUM_CORES=$1

# Get the output from the second argument
OUTPUT=$2

# Activate venv
eval "$(conda shell.bash hook)"

# Load new gwas_list
conda activate dashboard_project_env
echo "Loading files to process..."
python pipeline/gwas_list_manager.py
status=$?

# Check the exit status of the Python script
if [ $status -eq 1 ]; then
  exit 0
fi

# GenoPred dry run
cd ${genopred_dir}/pipeline
conda activate genopred

## Xavier TODO:I addedd this line because for some reason it went to look to the cache in the home directory and this can't happen in the cluster 
export XDG_CACHE_HOME=/imppc/labs/dnalab/xfarrer/.snakemake

echo "Starting GenoPred dry run..."
snakemake -n --configfile=input_gcat/config.yaml --use-conda $OUTPUT --unlock

if echo "$output" | grep -q "Nothing to be done"; then
  echo "Everything is up to date, no new GWAS detected."
  exit 0
fi

echo "Done."

cd ${polygenie_dir}

# Check headers
conda activate dashboard_project_env
echo "Checking headers..."
python pipeline/gwas_validator.py
exit_code=$?

# Check the exit code
if [ $exit_code -eq 1 ]; then
  exit 0
fi
echo "Done."

# Launch GenoPred
conda activate genopred
cd ${genopred_dir}/pipeline
snakemake -j${NUM_CORES} --configfile=input_gcat/config.yaml --use-conda $OUTPUT
#snakemake --profile /imppc/labs/dnalab/share/GEPETO/PolyGENIE_data/6-cluster_example --configfile=input_gcat/config.yaml --use-conda $OUTPUT

# Activate processing venv
cd ${polygenie_dir}
conda activate dashboard_project_env

# Copy profiles files
SOURCE_DIR="${genopred_dir}/pipeline/results_gcat/output/GCATcore/pgs/EUR/megaprs"
DEST_DIR="profiles_data"
SCRIPT="pipeline/main.py"
LOG_FILE="pipeline/logs/update_db_output.log"

echo "Saving new GenoPred files to PolyGenie/profiles_data..."
for dir in "$SOURCE_DIR"/*; do
    if [ -d "$dir" ]; then
        # Find and copy .profiles files to the destination directory
        find "$dir" -maxdepth 1 -type f -name "*.profiles" -exec cp {} "$DEST_DIR" \;
    fi
done

# Update db
#echo "Calculating regressions and prevalences..."
#python "$SCRIPT" 0 0 1 1 &> "$LOG_FILE" 
#echo "Pipeline completed!"

exit 0
