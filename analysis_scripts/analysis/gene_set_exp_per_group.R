#!/usr/bin/env Rscript

# Author: Lindsey Sydnor

# Purpose:
#   Take a seurat object and explore expression of gene sets of interest between
#   pre-identified groupings. Explored here are genes associated with
#   cerebellar malformations AND genes found on chromosome 21 PER germ layer

# TODO: make as script taking in path or turn into an importable utility.
#       Address TODOs within script.

my_packages <- c("Seurat", "Matrix", "glue", "ggplot2", "integration",
                 "patchwork", "escape", "pheatmap", "biomaRt")

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
        BiocManager::install(pkg, lib = .libPaths()[1], ask = FALSE)
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
csv_dir <- file.path(base_dir, "embryoid_output", "CSVs")

# load dataset
embryoid <- readRDS(file.path(obj_dir, "embryoid_germ_class.rds"))

### load gene sets to test

# get genes associated with cerebellar malformations
malf_genes <- read.csv(file.path(base_dir, "data", "references",
                                 "cbl_malf_genes.csv"))[[1]]

# use biomaRt to get human genes from ensembl database
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genes <- getBM(attributes = c("external_gene_name", "chromosome_name"),
               filters = "external_gene_name",
               values = rownames(embryoid),
               mart = ensembl)

# get list of chromosome 21 genes
c21_genes <- subset(genes, subset = (chromosome_name == 21))$external_gene_name

# initialize
gene_sets <- list("malf_genes" = malf_genes, "c21_genes" = c21_genes)

# gene list must be names list of gene sets to investigate
# obj must be seurat object
test_gs <- function(obj, gene_lists, groupings = "gs_1_germ_layer",
                    out_dir = file.path(image_dir, "enrichment"),
                    gene_fontsize = 3, assay = "SCT") {

  # extract names of each gene set
  gs_names <- names(gene_lists)

  # get original object name
  obj_name <- deparse(substitute(obj))

  i <- 0

  # for each gene set
  for (i in seq_along(gene_lists)) {
    print(i)
    # get name of gene set
    gs_name <- gs_names[i]
    
    # print message
    # print(glue("Running analysis on seurat obj",
    #            " '{deparse(substitute(obj))}' w/ gene set '{gs_name}'"))
    print(glue("Running analysis on seurat obj '{obj_name}' w/ gene set",
               " '{gs_name}'"))

    # create outdir for this gene set
    d <- file.path(out_dir, gs_name, obj_name)
    dir.create(d, showWarnings = FALSE)

    # find which genes from gene set are found in obj
    # Calculate the intersect
    intersected_genes <- intersect(rownames(obj), gene_lists[[i]])

    # Create the list with the named element
    included_genes <- list()
    included_genes[[gs_name]] <- intersected_genes

    # write those genes as outfile
    write.csv(included_genes,
              file.path(csv_dir,
              glue("{gs_name}_genes_included_in_{obj_name}.csv")),
              row.names = FALSE)

    # run enrichment analysis on these genes using ssGSEA method
    es <- enrichIt(obj = obj, gene.sets = included_genes, groups = 1000,
                   method = "ssGSEA")

    # add results as metadata
    obj[[glue("{gs_name}_enrich_ssGSEA")]] <- es

    # TODO why we also need this?
    obj <- AddMetaData(obj, es)

    # set as active ident
    Idents(obj) <- glue("{gs_name}_enrich_ssGSEA")

    obj@meta.data$active.idents <- obj@active.ident

    # create a metadata data frame including enrichment results
    es2 <- data.frame(obj[[]], Idents(obj))

    # create ridgeplot of gene set enrichment per groupings category
    p <- ridgeEnrichment(es2, gene.set = gs_name, group = groupings,
                         add.rug = TRUE)

    ggsave(file.path(d, glue("{groupings}_ridgeplot_ssGSEA.png")), plot = p)

    ### Create heatmaps of included genes (not gene set) per groupings
    # average expression for psuedobulking by germ_layer
    avg_exp <- AverageExpression(obj, features = included_genes[[gs_name]],
                                 return.seurat = FALSE, assays = assay,
                                 group.by = groupings, verbose = FALSE)[[1]]

    # create heatmap of pseudobulked expression in these genes
    p <- pheatmap(t(avg_exp), cluster_cols = FALSE,
                  fontsize_col = gene_fontsize)

    ggsave(plot = p, filename = file.path(d, "average_exp_heatmap.png"))

    ### find genes not expressed in specific germ layers

    # find array indices of 0 average expression (rows and column of avg_exp)
    arr_locs <- which(avg_exp == 0, arr.ind = TRUE)

    # if there is zero average expression for any of the genes of interest
    if (length(arr_locs) > 0) {
      # get names of columns (germ layer names) corresponding to array indices
      named_cols <- colnames(avg_exp[, which(avg_exp == 0, arr.ind = TRUE,
                                            useNames = TRUE)[, 2]])
      # convert to data frame
      arr_locs <- as.data.frame(arr_locs)
      # set map of data frame to not-expressed gene names and corresponding
      # layers
      col_mapping <- setNames(named_cols, rownames(arr_locs))
      # set germ layer column
      arr_locs$col <- col_mapping[rownames(arr_locs)]
      # extract only gene names and germ layer mapping (row number inessential)
      output <- arr_locs[2]

      # save output
      write.table(output, file.path(csv_dir,
                                    glue("{gs_name}_not_exp_w_layer_in_",
                                         "{obj_name}.csv")))
    }

    # TODO: append to info.txt, if already made in sequence, w/ description


    # make heatmap of comparative logFC per germ layer expression
    group_names <- unique(obj[[groupings]][[1]])

    # set idents to germ layer ID
    Idents(obj) <- groupings

    # initialize
    group_subsets <- list()

    # create cell subsets containing only germ layer of interest
    for (g in group_names) {
      group_subsets[[g]] <- subset(obj, ident = g)
    }

    # initialize enrichment matrix with germ layer names and genes of interest
    enrichment_matrix <- matrix(NA, nrow = length(group_names),
                                ncol = length(included_genes[[gs_name]]))
    rownames(enrichment_matrix) <- group_names
    colnames(enrichment_matrix) <- included_genes[[gs_name]]

    for (n in seq_along(group_names)) {
      # subset layer of interest
      group_subset <- group_subsets[[group_names[n]]]
      # iterate over genes of interest
      for (j in seq_along(included_genes[[gs_name]])) {
        gene <- included_genes[[gs_name]][j]
        # get average expression in this germ layer for this gene
        gene_expression <- mean(group_subset@assays[[assay]]@data[gene, ])
        # compare expression in this germ layer against average across all germ
        # layers
        enrichment_matrix[n, j] <- gene_expression / mean(avg_exp[gene, ])
      }
    }

    # some NANs were introduced, as there was no average expression for some
    # genes in specific germ layers
    p <- pheatmap(enrichment_matrix, cluster_rows = FALSE,
                  fontsize_col = gene_fontsize)

    # save
    ggsave(plot = p, filename = file.path(d, "FC_exp_heatmap.png"))
    # browser()
  }
}

# run on the gene sets of interest
test_gs(embryoid, gene_sets)

# TODO: change path if keeping
crouch <- readRDS("~/Downloads/crouch.rds")

test_gs(crouch, gene_sets, groupings = "cell_type", assay = "RNA",
        gene_fontsize = 7)
