#!/bin/bash
# request Bourne Again shell as shell for job
#$ -S /bin/bash
# Name for the script in the queuing system
#$ -N PolyGenie
# name of the queue you want to use
#$ -q d12imppc
#$ -cwd                             # Use the current working directory
#$ -o logs/process_profiles.log      # Standard output log file
#$ -e logs/process_profiles.err       # Standard error log file
#$ -t 1-<num_files>                 # Job array index, to be replaced later

# Load the Conda environment
# Replace 'myenv' with the name of your Conda environment
source /imppc/lab/dnalab/xfarrer/miniforge3/etc/profile.d/conda.sh  # Adjust this path to your Conda installation
conda activate dashboard_project_env  # Replace 'myenv' with the name of your Conda environment

# Create logs directory if it doesn't exist
mkdir -p logs

# Get the list of profile files
profiles_dir="profiles_data"
files=($(ls $profiles_dir/*.profiles))

# Get the number of files
num_files=${#files[@]}

# Generate the job array range for submission
if [ "$num_files" -gt 0 ]; then
    echo "Submitting job array for $num_files profile files."
else
    echo "No profile files found in $profiles_dir."
    exit 1
fi

# Prepare the command to be submitted as job array
for ((i=0; i<num_files; i++)); do
    file=${files[$i]}
    gwas_name=$(basename "$file" .profiles)  # Get the base name of the file without extension

    # Set the 'ori' value (you can customize this based on your logic)
    ori=0  # or set to 1 or any other value you need

    # Print the command to verify before submission
    echo "Submitting job for file: $file"
    echo "qsub -N process_${gwas_name} -b y -t 1-${num_files} python process_profile.py \"$file\" \"$ori\""
done

# Now, submit the job array command with the calculated range
qsub -N process_profiles -t 1-$num_files -b y python process_profile.py "${files[@]}" "$ori"

