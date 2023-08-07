#!/usr/bin/env Rscript

# Purpose:
#   TODO

# THIS DOESN'T WORK
Sys.setenv(R_INTERACTIVE = "false")

# load packages
library(Seurat)
library(dplyr)
library(patchwork)
library(glue)
library(ggplot2)
library(clustree)
library(ggraph)
library(pryr)
library(integration) # made in utilities

############################ Record resource usage #############################

# Record the starting memory usage
start_mem <- mem_used()
# Record the starting time
start_time <- proc.time()

###################### Set up background (dirs, presets) ######################

# RELATIVE PATH for github (call script anywhere within repo file structure):
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

############################## QC investigation ##############################

# set active assay to RNA
DefaultAssay(organoid) <- "RNA"

# set out dir
d <- file.path(umap_dir, "pre_QC")
dir.create(d, showWarnings = FALSE)

# plot original UMAP
p <- DimPlot(organoid, reduction = "umap") + NoAxes()
ggsave(filename = file.path(d, "cell_classes.png"), plot = p)

# explore metadata
p <- DimPlot(organoid, reduction = "umap", group.by = "stim") + NoAxes()
ggsave(filename = file.path(d, "stim.png"), plot = p)

# set out dir
d <- file.path(image_dir, "QC", "pre_rem", "violin_plots")
dir.create(d, recursive = TRUE, showWarnings = FALSE)

# create violin plots of features of interest per cell type
to_test <- c("percent.mt", "nUMI", "nFeature_RNA", "nCount_RNA")
for (test in to_test) {
  p <- VlnPlot(organoid, features = test)
  ggsave(plot = p, filename = file.path(d, glue("{test}.png")))
}

# dataset came excluding cells with < 1285 unique genes... odd
# min(organoid$nFeature_RNA) # 1285

# create elbowplot
p <- ElbowPlot(organoid)
ggsave(filename = file.path(d, "elbow_plot.png"), plot = p)

################################ Preprocessing ################################
# retain cell class names in metadata
cell_classes <- Idents(organoid)
organoid$cell_classes <- gsub(" ", "_", cell_classes) # rid names of spaces

# combined 6 organoid datasets, so run on "integrated"
DefaultAssay(organoid) <- "integrated"

if (length(organoid@graphs) == 0) {
  print("No precomputed graphs in this object. Creating NN and SNN graphs...")
  organoid <- FindNeighbors(organoid, reduction = "pca", dims = 1:npcs)
}

# rerun processing pipeline @ multiple resolutions on integrated assay:
resolutions <- seq(0.1, 1, by = 0.1)

# set outdir
d <- file.path(umap_dir, "pre_QC", "by_res")
dir.create(d, showWarnings = FALSE)

# Find clusters @ each resolution
for (r in resolutions) {
  organoid <- FindClusters(organoid, verbose = FALSE, resolution = r)
  p <- DimPlot(organoid, reduction = "umap",
               group.by = glue("integrated_snn_res.{r}")) + NoAxes()
  ggsave(filename = file.path(d, glue("{r}.png")), plot = p)
}

# save checkpoint
saveRDS(object = organoid, file = file.path(obj_dir, "organoid_updated.rds"))

# run clustree for ideal resolution info
clust_diagram <- clustree(organoid, prefix = "integrated_snn_res.")
ggsave(filename = file.path(image_dir, "clustree.png"), plot = clust_diagram)

############ Cell removal -- ssve RDS for optional future analysis ############

# thresholding -- excluding cells > 10% mito RNA
organoid_rem <- subset(x = organoid, subset = percent.mt <= 10)

### remake QC / UMAP plots
# set out dir
d <- file.path(umap_dir, "post_QC")
dir.create(d, showWarnings = FALSE)

# plot original UMAP
p <- DimPlot(organoid_rem, reduction = "umap") + NoAxes()
ggsave(filename = file.path(d, "cell_classes.png"), plot = p)

# explore metadata
p <- DimPlot(organoid_rem, reduction = "umap", group.by = "stim") + NoAxes()
ggsave(filename = file.path(d, "stim.png"), plot = p)

# set out dir
d <- file.path(image_dir, "QC", "post_rem", "violin_plots")
dir.create(d, recursive = TRUE, showWarnings = FALSE)

# create violin plots of features of interest per cell type
to_test <- c("percent.mt", "nUMI", "nFeature_RNA", "nCount_RNA")
for (test in to_test) {
  p <- VlnPlot(organoid_rem, features = test)
  ggsave(plot = p, filename = file.path(d, glue("{test}.png")))
}

saveRDS(organoid_rem, file.path(obj_dir, "organoid_QC.rds"))

# write info file to explain organization of obj_dir (overwrites currently)
file_name <- "info.txt"
line <-
  "organoid_QC: \
  \tNot yet used in analysis other than to create images 'pre_QC' or 'pre_rem'."
# Write the data to the file
writeLines(line, file.path(obj_dir, file_name)) # append = TRUE if not overwrite


############################ Output resource usage #############################

# Record the ending time
end_time <- proc.time()
# Calculate the elapsed time
elapsed_time <- end_time - start_time

# Record the ending memory usage
end_mem <- mem_used()
# Calculate the memory usage during execution
memory_used <- end_mem - start_mem

# Print the computed resources
print(paste("Elapsed Time: ", format(elapsed_time[3], digits = 2), " seconds"))
print(paste("Memory Used: ", format(memory_used / 1024^2, digits = 2), " MB"))