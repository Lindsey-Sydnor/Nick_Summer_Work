#!/usr/bin/env Rscript

# Author: Lindsey Sydnor

# Purpose:
#   Take a seurat object and explore expression of cerebellar malformation genes

# TODO: make as script taking in path or turn into an importable utility.

my_packages <- c("Seurat", "Matrix", "glue", "ggplot2", "integration",
                 "patchwork", "escape", "pheatmap")

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

# get genes associated with cerebellar malformations
malf_genes <- read.csv(file.path(base_dir, "data", "references",
                                 "cbl_malf_genes.csv"))[[1]]

# find which cerebellar malformation genes are found in embryoid dataset
included_genes <- list("malf_genes" = intersect(rownames(embryoid), malf_genes))

# save
write.csv(included_genes,
          file.path(csv_dir, "cblr_malf_genes_included.csv"),
          row.names = FALSE)

# run enrichment analysis on these genes using ssGSEA method
es <- enrichIt(obj = embryoid, gene.sets = included_genes, groups = 1000,
               method = "ssGSEA")

# check not all 20,000 genes (check that it's specific to dataset composition)

# add enrichment metadata
embryoid$malf_enrich_ssGSEA <- es
embryoid <- AddMetaData(embryoid, es)
# Idents(embryoid) <- "malf_enrich_ssGSEA"
embryoid@meta.data$active.idents <- embryoid@active.ident

es2 <- data.frame(embryoid[[]], Idents(embryoid))

# create ridgeplot of enrichment per germ layer category
p <- ridgeEnrichment(es2, gene.set = "malf_genes", group = "gs_1_germ_layer",
                     add.rug = TRUE)

d <- file.path(image_dir, "enrichment")
dir.create(d, showWarnings = FALSE)

ggsave(file.path(d, "GS1_ridgeplot_ssGSEA.png"), plot = p)

# run enrichment analysis on these genes using UCell method
es <- enrichIt(obj = embryoid, gene.sets = included_genes, groups = 1000,
               method = "UCell")

# add enrichment metadata
embryoid$malf_enrich_UCell <- es
Idents(embryoid) <- "malf_enrich_UCell"

# prep input
es2 <- data.frame(embryoid[[]], Idents(embryoid))

# create ridgeplot of enrichment per germ layer category
p <- ridgeEnrichment(es2, gene.set = "malf_genes", group = "gs_1_germ_layer",
                     add.rug = TRUE)

ggsave(file.path(d, "GS1_ridgeplot_UCell.png"), plot = p)




# heatmaps per gene (not gene set)

# average expression for psuedobulking by germ_layer
avg_exp <- AverageExpression(embryoid, features = included_genes[[1]],
                             return.seurat = FALSE, assays = "SCT",
                             group.by = "gs_1_germ_layer", verbose = FALSE)[[1]]


# create heatmap of pseudobulked expression in these genes
p <- pheatmap(avg_exp[[1]], cluster_rows = FALSE, fontsize = 8)

ggsave(plot = p, filename = file.path(d, "average_exp_heatmap.png"))


arr_locs <- which(avg_exp == 0, arr.ind = TRUE)
named_cols <- colnames(avg_exp[,which(avg_exp == 0, arr.ind = TRUE,
                               useNames = TRUE)[,2]])
arr_locs <- as.data.frame(arr_locs)
col_mapping <- setNames(named_cols, rownames(arr_locs))
arr_locs$col <- col_mapping[rownames(arr_locs)]

output <- arr_locs[2]

print(output)

# check not all 20,000 genes (check that it's specific to dataset composition)