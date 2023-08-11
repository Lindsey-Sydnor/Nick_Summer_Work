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

gs1 <- list("Endoderm" = c("FOXA2", "SOX17", "GATA6", "AFP"),
            "Mesoderm" = c("MESP1", "EOMES", "TBXT", "MIXL1", "HAND1",
                           "HAND2", "TBX6"),
            "Ectoderm" = c("NEUROD1", "PAX6", "SOX1", "SALL3"),
            "Pluripotent" = c("POU5F1", "SOX2"))

# load in embryoid object
embryoid <- readRDS(file.path(obj_dir, "embryoid_update.rds"))

# create visualizations for each gene set
gene_set_visualizations <- function(seurat_obj, gene_sets, gs_names = NULL) {
  i <- 1 # initialize
  for (gs in gene_sets) {
    if (is.null(gs_names)) {
      # if gene sets didn't come with associated names, use numerical naming
      # scheme
      gs_name <- i
    } else {
      # if gene sets come with associated names, name the outfiles with them
      gs_name <- gs_names[i]
    }
    # create outdir for germ-layer-specific analysis
    germ_dir <- file.path(image_dir, "germ_layers",
                          "germ_layer_marker_exploration",
                          glue("gene_set_{gs_name}"))
    dir.create(germ_dir, recursive = TRUE, showWarnings = FALSE)

    # format gene set for readability to write to info.txt in germ_dir
    formatted_text <- paste(
      "Exploring feasibility of gene set:\n\n",
      paste(sapply(names(gs), function(name) {
        paste0("$", name, "\n[1] \"", paste(gs[[name]], collapse = "\" \""),
                "\"")
        }), collapse = "\n\n")
    )

    # Write the text to the file
    file_name <- "info.txt"
    writeLines(formatted_text, file.path(germ_dir, file_name))

    # initialize
    included_GMs <- list()

    # for each named germ_layer within gene set of interest
    for (germ_layer in names(gs)) {
      # determine which genes are present in seurat_obj
      matching_genes <- intersect(rownames(seurat_obj),
                                  gs[[germ_layer]])
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
      seurat_obj <- AddModuleScore(seurat_obj,
                                   name = glue("gs_{gs_name}_{germ_layer}"),
                                   features = list(included_GMs[[germ_layer]]),
                                   assay = "SCT")

      # create FeaturePlot of module score for this germ layer's genes
      p <- FeaturePlot(seurat_obj,
                       features = glue("gs_{gs_name}_{germ_layer}1"),
                       repel = TRUE) +
        labs(caption = glue("Genes used to define set: ",
                            "{list(included_GMs[[germ_layer]])}")) +
        theme(plot.caption = element_text(hjust = 0.5)) # center caption

      # save plot
      ggsave(filename = file.path(d, glue("{germ_layer}_unicolor.png")),
             plot = p)

      # repeat, but with two tones
      p <- p +
        scale_colour_gradientn(colours = rev(brewer.pal(n = 4, name = "RdBu")))

      # save plot
      ggsave(filename = file.path(d, glue("{germ_layer}_multicolor.png")),
             plot = p)
    }

    # Define a function to assign categories
    assign_category <- function(row_values) {
      category <- rep("unknown", nrow(row_values))

      pluri_val <- row_values[glue("gs_{gs_name}_Pluripotent1")]
      endo_val <- row_values[glue("gs_{gs_name}_Endoderm1")]
      ecto_val <- row_values[glue("gs_{gs_name}_Ectoderm1")]
      meso_val <- row_values[glue("gs_{gs_name}_Mesoderm1")]

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
    metadata_cols <- c(glue("gs_{gs_name}_Pluripotent1"),
                       glue("gs_{gs_name}_Endoderm1"),
                       glue("gs_{gs_name}_Ectoderm1"),
                       glue("gs_{gs_name}_Mesoderm1"))
    values_matrix <- seurat_obj@meta.data[metadata_cols]

    # Apply the function to each row (cell) of the values matrix
    categories <- assign_category(values_matrix)

    # add a metadata column to seurat_obj containing this classification
    seurat_obj[[glue("gs_{gs_name}_germ_layer")]] <- categories

    p <- DimPlot(seurat_obj, group.by = glue("gs_{gs_name}_germ_layer"),
                 reduction = "umap", pt.size = 1) +
      ggtitle(glue("Germ Layer Classificatio via Gene Set {gs_name}")) +
      NoAxes()

    ggsave(filename = file.path(germ_dir, "germ_layer_class.png"), plot = p)

    # save counts of each germ layer to csv
    write.csv(table(seurat_obj[[glue("gs_{gs_name}_germ_layer")]]),
              file = file.path(germ_dir, "germ_layer_class.csv"),
              row.names = FALSE)

    i <- i + 1 # increment counter
  }
  return(seurat_obj)
}

# create a function to combine gene sets and get unique genes for each layer
combine_and_unique <- function(genes1, genes2) {
  unique(c(genes1, genes2))
}

# get a combined list of unique genes in gene set of interest + Kim's gene
# set using vector operations
all_GMs <- Map(combine_and_unique, kim_GMs, gs1)

# call fxn
gs_vec <- list(gs1, all_GMs)
gs_names <- list(1, "1_and_Kim")
embryoid <- gene_set_visualizations(embryoid, gs_vec, gs_names)

### TODO: something about covarying genes? DE expression for identified cell
# classes?