#!/usr/bin/env Rscript

# Purpose:
#   TODO

################# Set up background (dirs, presets) #################

# load packages
library(Seurat)
library(dplyr)
library(patchwork)
library(glue)
library(ggplot2)
library(clustree)
library(ggraph)

# RELATIVE PATH for github (call script anywhere within repo file structure):
initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name,
                   initial.options)])
# Get the directory of the script
script.dir <- dirname(normalizePath(script.name))
# Get the parent directory of the script (one level up)
parent.dir <- dirname(script.dir)
# Get the parent directory of the parent directory (two levels up)
base_dir <- dirname(parent.dir)

# creating output directories
dirs <- c("images", "objs")
for (d in dirs){
  dir.create(file.path(base_dir, d), recursive = TRUE, showWarnings = FALSE)
}

# setting defaults
image_dir <- file.path(base_dir, "images")
obj_dir <- file.path(base_dir, "objs")
data_dir <- file.path(base_dir, "data")
umap_dir <- file.path(image_dir, "umaps")

# create dirs, if not already made
dir.create(umap_dir, showWarnings = FALSE)

# presets
npcs <- 30

# load object
organoid <- readRDS(file.path(data_dir, "organoid.rds"))

# NOTE: PCA, UMAP reductions precomputed. Small dataset (1653 cells).

######################## QC + investigation ########################

# set active assay to RNA
DefaultAssay(organoid) <- "RNA"

# plot original UMAP
p <- DimPlot(organoid, reduction = "umap") + NoAxes()
ggsave(filename = file.path(umap_dir, "orig.png"), plot = p)

# explore metadata
p <- DimPlot(organoid, reduction = "umap", group.by = "stim") + NoAxes()
ggsave(filename = file.path(umap_dir, "stim.png"), plot = p)


# create QC and vln plot dir
d <- file.path(image_dir, "QC", "violin_plots")
dir.create(d, recursive = TRUE, showWarnings = FALSE)

# create violin plots of features of interest per cell type
to_test <- c("percent.mt", "nUMI", "nFeature_RNA", "nCount_RNA")
for (test in to_test) {
  p <- VlnPlot(organoid, features = test)
  ggsave(plot = p, filename = file.path(d, glue("{test}.png")))
}

# dataset came excluding cells with < 1285 unique genes... odd
min_feat <- min(organoid$nFeature_RNA) # 1285

######################## Reconstruct Seurat object ########################

# thresholding -- excluding cells > 5% mito RNA
organoid <- subset(x = organoid, subset = percent.mt <= 5)

# perform preprocesing
organoid <- FindVariableFeatures(organoid)

# retain cell class names in metadata (not currently included)
cell_classes <- Idents(organoid)
organoid$cell_classes <- gsub(" ", "_", cell_classes) # rid names of spaces

# create elbowplot
p <- ElbowPlot(organoid)
ggsave(filename = file.path(image_dir, "elbow_plot.png"), plot = p)

# rerun processing pipeline @ multiple resolutions on RNA assay:
resolutions <- seq(0.1, 1, by = 0.1)

if (length(organoid@graphs) == 0) {
  print("No precomputed graphs in this object. Creating NN and SNN graphs...")
}

# create graphs
organoid <- FindNeighbors(organoid, dims = 1:npcs, assay = "RNA",
                          verbose = FALSE)