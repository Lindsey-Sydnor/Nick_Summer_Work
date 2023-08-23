#!/bin/bash

# Get the directory of the currently executing script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Construct script paths relative to this script's directory
R_SCRIPTS_PATH="${SCRIPT_DIR}/analysis_scripts/analysis/"

# make sure all of the scripts are executable
chmod +x $R_SCRIPTS_PATH/*.R

# load R module
module load R/4.1.0-foss-2020b

# Call scripts in the desired sequence
echo "Installing utilities..."
"${R_SCRIPTS_PATH}create_package.R"

echo "Retrieving data..."
"${SCRIPT_DIR}/data/data_accession.sh"

echo "Running paper's analysis..."
qsub -P $(Project code AutismR21) -q paidq -l walltime=3:00:00 -l \
  select=1:ncpus=1:mem=64gb \
  "$R_SCRIPTS_PATH/PBS_scripts/run_paper_analysis.pbs"

# FLUSH OUT EVENTUALLY
# echo "Processing data..."
# "${R_SCRIPTS_PATH}data_processing.R"

# echo "Running analysis..."
# "${R_SCRIPTS_PATH}analysis.R"

# echo "Generating plots..."
# "${R_SCRIPTS_PATH}plotting.R"

echo "Analysis pipeline completed!"
