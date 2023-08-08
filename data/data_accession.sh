#!/bin/bash

# Setup directories

# Determine the directory of the script
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"

# Define data outdirs
OUTS=("embryoid" "organoid")

# Make data outdirs
for D in "${OUTS[@]}"; do
  mkdir "$SCRIPT_DIR/$D"
done


# Download GEO data for organoid

# Set BASE_DIR to the organoid out directory
BASE_DIR="$SCRIPT_DIR/organoid"

# Define URL
URL="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150153/suppl/\
GSE150153_RAW.tar"

# Use wget to download the dataset and save it to the specified folder
wget -P "$BASE_DIR" "$URL"

# Untar and gunzip everything in organoid dir
tar -xvf $BASE_DIR/* -C $BASE_DIR

# Remove the tar
rm $BASE_DIR/GSE150153_RAW.tar

# Download GEO data for embryoid

# Set BASE_DIR to the embryoid out directory
BASE_DIR="$SCRIPT_DIR/embryoid"

# Download control-specific URLs
URLS=(
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3573nnn/GSM3573649/suppl/\
GSM3573649_D_barcodes.tsv.gz"
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3573nnn/GSM3573649/suppl/\
GSM3573649_D_filtered_gene_bc_matrices_h5.h5"
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3573nnn/GSM3573649/suppl/\
GSM3573649_D_genes.tsv.gz"
    "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3573nnn/GSM3573649/suppl/\
GSM3573649_D_matrix.mtx.gz"
)

# Loop through the URLs and download datasets to the specified folder
for URL in "${URLS[@]}"; do
    wget -P "$BASE_DIR" "$URL"
done

# gunzip everything in embryoid dir
gunzip $BASE_DIR/*.gz
