#!/usr/bin/env Rscript

# Author: Lindsey Sydnor

# Purpose:
#   Determine thresholds used for preprocessing based on raw vs. processed
#   object.

################# Set up background (dirs, presets) #################

# get location of base_dir irrespective of where script is called
initial_options <- commandArgs(trailingOnly = FALSE)
file_arg_name <- "--file="
script_name <- sub(file_arg_name, "",
                   initial_options[grep(file_arg_name,
                                        initial_options)])
# Get the directory of the script
script_dir <- dirname(normalizePath(script_name))
# Get the parent directory of the script (one level up)
parent_dir <- dirname(script_dir)
# Get the parent directory of the parent directory (two levels up)
base_dir <- dirname(parent_dir)

# setting defaults
data_dir <- file.path(base_dir, "data")

# load in objects
updated_organoid <- readRDS(file.path(data_dir, "updated_organoid.rds"))
raw_organoid <- readRDS(file.path(data_dir, "organoid.rds"))

# determine max percent.mt
raw_organoid@active.assay <- "RNA" # look @ same assay
max(raw_organoid$percent.mt) # 38.34576
max(updated_organoid$percent.mt) # 38.34576

# Didn't differ. High percent.mt cells were not excluded from the RDS's we have
# available. Means we need to reconstruct preprocessed RDS, set new thresholds
# that may differ from past analyses. Need to reconstruct figures.