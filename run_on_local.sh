#!/bin/bash

# Get the directory of the currently executing script
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# Construct script paths relative to this script's directory
R_SCRIPTS_PATH="${SCRIPT_DIR}/analysis_scripts/analysis/"

# Call scripts in the desired sequence
echo "Installing utilities..."
"${R_SCRIPTS_PATH}create_package.R"

echo "Retrieving data..."
"${SCRIPT_DIR}/data/data_accession.sh"

# FLUSH OUT EVENTUALLY

# echo "Processing data..."
# "${R_SCRIPTS_PATH}data_processing.R"

# echo "Running analysis..."
# "${R_SCRIPTS_PATH}analysis.R"

# echo "Generating plots..."
# "${R_SCRIPTS_PATH}plotting.R"

echo "Analysis pipeline completed!"
