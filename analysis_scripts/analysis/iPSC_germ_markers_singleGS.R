#!/usr/bin/env Rscript

# Author: Lindsey Sydnor

# Purpose:
#   An exploration of potential pluripotent and germ layer marker genes for the
#   embryoid dataset.

my_packages <- c("Seurat", "glue", "ggplot2", "integration", "patchwork",
                 "RColorBrewer")

# define repo to install from
options(repos = c(CRAN = "http://cran.rstudio.com"))

# install / import using CRAN or Bioconductor, as appropriate
for (pkg in my_packages) {
  if (!require(pkg, character.only = TRUE)) {
    tryCatch(
      {install.packages(pkg, lib = .libPaths()[1])
       invisible(lapply(pkg, library, character.only = TRUE))},
      error = function(e) {
        message("Installation via install.packages() failed. Trying
                 BiocManager.")
        BiocManager::install(pkg, lib = .libPaths()[1])
        invisible(lapply(pkg, library, character.only = TRUE))
      },
      warning = function(w) {
        message(w)
      }
    )
  }
}

base_dir <- getwd() # TODO: CHANGE FOR GIT
obj_dir <- file.path(base_dir, "embryoid_output", "objs")
image_dir <- file.path(base_dir, "embryoid_output", "images")

# potential germ layer gene markers
kim_GMs <- list("Endoderm" = c("SOX17", "FOXA2", "CXCR4", "GATA4"),
                "Mesoderm" = c("NCAM1", "TBXT"),
                "Ectoderm" = c("NES", "PAX6"),
                "Pluripotent" = c())

pruned_GMs <- list("Endoderm" = c("FOXA2", "SOX17", "GATA6", "AFP"),
                   "Mesoderm" = c("MESP1", "EOMES", "TBXT", "MIXL1", "HAND1",
                                  "HAND2", "TBX6"),
                   "Ectoderm" = c("NEUROD1", "PAX6", "SOX1", "SALL3"),
                   "Pluripotent" = c("POU5F1", "SOX2"))

### use vector operations to produce master list
# create a function to combine and unique genes for each germ layer
combine_and_unique <- function(genes1, genes2) {
  unique(c(genes1, genes2))
}

# use Map to apply the function to each germ layer
all_GMs <- Map(combine_and_unique, kim_GMs, pruned_GMs)

# create outdir for germ-layer-specific analysis
germ_dir <- file.path(image_dir, "germ_layers", "germ_layer_marker_exploration")
dir.create(germ_dir, recursive = TRUE, showWarnings = FALSE)

# load in embryoid object
embryoid <- readRDS(file.path(obj_dir, "embryoid_update.rds"))

# initialize
included_GMs <- list()

for (germ_layer in names(pruned_GMs)) {
  # determine which genes are present in embryoids
  matching_genes <- intersect(rownames(embryoid), pruned_GMs[[germ_layer]])
  # store genes present in dataset
  included_GMs[[germ_layer]] <- matching_genes
}

# initialize
all_markers <- list() # to store the markers for each germ layer
threshold <- 0 # set threshold for expression (0)

# Define your germ layer gene sets
germ_layers <- names(included_GMs)

# Loop through each germ layer
for (germ_layer in germ_layers) {
  # make an outdir to save to for each layer
  d <- file.path(germ_dir, "feature_plots")
  dir.create(d, showWarnings = FALSE)

  # create a binary module score indicating if the cell expresses the germ
  # layer's genes
  embryoid <- AddModuleScore(embryoid, name = germ_layer,
                             features = list(included_GMs[[germ_layer]]),
                             assay = "SCT")

  # create FeaturePlot of module score for this germ layer's genes
  p <- FeaturePlot(embryoid, features = glue("{germ_layer}1"),
                   repel = TRUE) +
    labs(caption = glue("Genes used to define set: ",
                        "{list(included_GMs[[germ_layer]])}")) +
    theme(plot.caption = element_text(hjust = 0.5)) # center caption

  # save plot
  ggsave(filename = file.path(d, glue("{germ_layer}_unicolor.png")), plot = p)

  # repeat, but with two tones
  p <- p +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 4, name = "RdBu")))

  # save plot
  ggsave(filename = file.path(d, glue("{germ_layer}_multicolor.png")), plot = p)
}

# Define a function to assign categories
assign_category <- function(row_values) {  
  category <- rep("unknown", nrow(row_values))

  pluri_val <- row_values["Pluripotent1"]
  endo_val <- row_values["Endoderm1"]
  ecto_val <- row_values["Ectoderm1"]
  meso_val <- row_values["Mesoderm1"]

  condition <- rowSums(row_values > 0.2) >= 2
  category[condition] <- "conflicting"

  condition <- pluri_val > 0.2 & rowSums(row_values > 0.2) < 2
  category[condition] <- "pluripotent"

  condition <- endo_val > 0.2 & rowSums(row_values > 0.2) < 2
  category[condition] <- "endoderm"

  condition <- ecto_val > 0.2 & rowSums(row_values > 0.2) < 2
  category[condition] <- "ectoderm"

  condition <- meso_val > 0.2 & rowSums(row_values > 0.2) < 2
  category[condition] <- "mesoderm"

  return(category)
}

# Extract the relevant metadata columns
metadata_cols <- c("Pluripotent1", "Endoderm1", "Ectoderm1", "Mesoderm1")
values_matrix <- embryoid@meta.data[metadata_cols]

# Apply the function to each row (cell) of the values matrix
categories <- assign_category(values_matrix)

# add a metadata column to embryoid containing this classification
embryoid$germ_layer <- categories

p <- DimPlot(embryoid, group.by = "germ_layer", reduction = "umap", pt.size = 1) +
  ggtitle("Germ Layer Classification") + NoAxes() #+
  # scale_color_manual(values = brewer.pal(6, "Set2"))

ggsave(filename = file.path(base_dir, "test.png"), plot = p)


### TODO: something about covarying genes