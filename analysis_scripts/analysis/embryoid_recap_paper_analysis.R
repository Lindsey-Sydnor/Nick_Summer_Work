#!/usr/bin/env Rscript

# Author: Lindsey Sydnor

# Purpose:
#   Meant to be an exact replication of embryoid body dataset analysis as in
#   paper (not my chosen presets)

# TODO: uncomment sections performing preprocesing for final upload

my_packages <- c("Seurat", "Matrix", "glue", "ggplot2", "integration")

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

# presets
npcs <- 20

# set paths
# base_dir <- getwd() # TODO: CHANGE FOR GIT
base_dir <- "/active/aldinger_k/kimslab_temp/test_nick/git_repo/Nick_Summer_Work"
data_dir <- file.path(base_dir, "data", "embryoid") # premade, stores raw data
obj_dir <- file.path(base_dir, "embryoid_output", "paper_analysis", "objs")
image_dir <- file.path(base_dir, "embryoid_output", "paper_analysis", "images")
csv_dir <- file.path(base_dir, "embryoid_output", "paper_analysis", "CSVs")

# creating output directories
dirs <- c(image_dir, obj_dir, csv_dir)
for (d in dirs){
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# # read in data for control
# counts_ctl <- readMM(file.path(data_dir, "GSM3573649_D_matrix.mtx"))
# # TODO: ^ large file, in .gitignore for now
# barcodes_ctl <- read.table(file.path(data_dir, "GSM3573649_D_barcodes.tsv"))
# genes_ctl <- read.table(file.path(data_dir, "GSM3573649_D_genes.tsv"))

# # construct counts matrix for Seurat
# colnames(counts_ctl) <- barcodes_ctl$V1
# rownames(counts_ctl) <- genes_ctl$V2 # extract gene symbols, ignore ensembl IDs

# # TODO: interesting... non-unique features
# embryoid_ctl <- CreateSeuratObject(counts = counts_ctl,
#                                    project = "embryoid_ctl",
#                                    min.features = 200)

# saveRDS(embryoid_ctl, file.path(data_dir, "embryoid_control.rds"))
embryoid_ctl <- readRDS(file.path(data_dir, "embryoid_control.rds"))

# read in data for nicotine
counts_nic <- readMM(file.path(data_dir, "GSM3573650_N_matrix.mtx"))
# TODO: ^ large file, in .gitignore for now
barcodes_nic <- read.table(file.path(data_dir, "GSM3573650_N_barcodes.tsv"))
genes_nic <- read.table(file.path(data_dir, "GSM3573650_N_genes.tsv"))

# construct counts matrix for Seurat
colnames(counts_nic) <- barcodes_nic$V1
rownames(counts_nic) <- genes_nic$V2 # extract gene symbols, ignore ensembl IDs

# TODO: interesting... non-unique features
embryoid_nic <- CreateSeuratObject(counts = counts_nic,
                                   project = "embryoid_nic",
                                   min.features = 200)

# saveRDS(embryoid_nic, file.path(data_dir, "embryoid_nicotine.rds"))
embryoid_nic <- readRDS(file.path(data_dir, "embryoid_nicotine.rds"))

# embryoid_ctl <- subset(embryoid_ctl, subset = nFeature_RNA < 6000)

# # subset out cells with high mitochondrial content
# embryoid_ctl <- PercentageFeatureSet(embryoid_ctl, pattern = "^MT-",
#                                      col.name = "percent.mt")
# embryoid_ctl <- subset(embryoid_ctl, subset = percent.mt < 20)
# ncol(embryoid_ctl) == 6766 # the num in their paper

# embryoid_ctl <- NormalizeData(embryoid_ctl,
#                               normalization.method = "LogNormalize",
#                               scale.factor = 10000)

# # Filter genes based on average expression, dispersion, and cutoffs
# min_avg_expr <- 0.0125
# max_avg_expr <- 3
# min_dispersion <- 0.5

# embryoid_ctl <- FindVariableFeatures(object = embryoid_ctl,
#                                      mean.function = ExpMean,
#                                      dispersion.function = LogVMR,
#                                      x.low.cutoff = min_avg_expr,
#                                      x.high.cutoff = max_avg_expr,
#                                      y.cutoff = min_dispersion)

# embryoid_ctl <- ScaleData(embryoid_ctl, features = rownames(embryoid_ctl),
#                           vars.to.regress = c("percent.mt", "nCount_RNA"))

# embryoid_ctl <- RunPCA(embryoid_ctl, assay = "RNA",
#                        features = VariableFeatures(embryoid_ctl),
#                        npcs = npcs)

# embryoid_ctl <- FindNeighbors(embryoid_ctl, dims = 1:npcs, assay = "RNA",
#                               verbose = FALSE)

# embryoid_ctl <- FindClusters(embryoid_ctl, verbose = FALSE, resolution = 0.8)

# # 20 dims
# embryoid_ctl <- RunTSNE(embryoid_ctl, dims = 1:npcs)

# p <- DimPlot(embryoid_ctl, reduction = "tsne") + NoAxes() +
#   ggtitle("tSNE Control")

# ggsave(file.path(image_dir, "umap.png"), plot = p)

# saveRDS(embryoid_ctl, file.path(obj_dir, "embryoid_control.rds"))


# markers_ctl <- FindMarkers(embryoid_ctl, logfc.threshold = 0.25, min.pct = 0.25,
#                            only.pos = TRUE)

# write.csv(markers_ctl, file.path(csv_dir, "FindMarkers_control.rds"))
# # then threshold for p val?

perform_seurat_analysis <- function(seurat_list, npcs, image_dir, obj_dir,
                                    csv_dir, res = 0.8) {
  for (i in seq_along(seurat_list)) {
    seurat_obj <- seurat_list[[i]]
    # # get original name of seurat object
    seurat_name <- names(seurat_list[i])

    # Subset cells based on nFeature_RNA
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA < 6000)
    
    # Remove cells with high mitochondrial content
    seurat_obj <- PercentageFeatureSet(seurat_obj, pattern = "^MT-",
                                       col.name = "percent.mt")
    seurat_obj <- subset(seurat_obj, subset = percent.mt < 20)
    
    # Normalize data
    seurat_obj <- NormalizeData(seurat_obj,
                                normalization.method = "LogNormalize",
                                scale.factor = 10000)
    
    # Filter genes based on average expression, dispersion, and cutoffs
    min_avg_expr <- 0.0125
    max_avg_expr <- 3
    min_dispersion <- 0.5
    
    # seurat_obj <- FindVariableFeatures(object = seurat_obj,
    #                                    mean.function = ExpMean,
    #                                    dispersion.function = LogVMR,
    #                                    x.low.cutoff = min_avg_expr,
    #                                    x.high.cutoff = max_avg_expr,
    #                                    y.cutoff = min_dispersion)

    seurat_obj <- FindVariableFeatures(object = seurat_obj,
                                    mean.function = ExpMean,
                                    dispersion.function = LogVMR,
                                    mean.cutoff = c(min_avg_expr, max_avg_expr),
                                    dispersion.cutoff = c(min_dispersion, Inf))
    
    # Scale data
    seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj),
                            vars.to.regress = c("percent.mt", "nCount_RNA"))
    
    # Run PCA
    seurat_obj <- RunPCA(seurat_obj, assay = "RNA",
                         features = VariableFeatures(seurat_obj), npcs = npcs)
    
    # Find neighbors
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:npcs, assay = "RNA",
                                verbose = FALSE)
    
    # Find clusters
    seurat_obj <- FindClusters(seurat_obj, verbose = FALSE, resolution = res)
    
    # Run tSNE
    seurat_obj <- RunTSNE(seurat_obj, dims = 1:npcs)
    
    # Plot tSNE plot
    p <- DimPlot(seurat_obj, reduction = "tsne") + NoAxes() +
      ggtitle(paste("tSNE ", seurat_name))

    # save
    ggsave(file.path(image_dir, paste0("tsne_", seurat_name, ".png")), plot = p)
    
    # Save Seurat object with seurat_name in the file name
    saveRDS(seurat_obj, file.path(obj_dir, paste0("embryoid_", seurat_name,
                                                  ".rds")))
    seurat_obj <- readRDS(file.path(obj_dir, paste0("embryoid_", seurat_name,
                                                    ".rds")))
    
    for (clust in unique(seurat_obj[[glue("RNA_snn_res.{res}")]])[[1]]) {
      # Find markers
      markers <- FindMarkers(seurat_obj, ident.1 = toString(clust),
                             logfc.threshold = 0.25, min.pct = 0.25,
                             only.pos = TRUE)

      d <- file.path(csv_dir, glue("{seurat_name}_FindMarkers"),
                     glue("res{res}"))
      dir.create(d, showWarnings = FALSE, recursive = TRUE)

      # save as CSV
      write.csv(markers, file.path(d, glue("clust{clust}.csv")))
    }
s
    pow_dir <- file.path(image_dir, "power_plots", seurat_name,
                         glue("res{res}"))
    dir.create(pow_dir, showWarnings = FALSE, recursive = TRUE)

    # make power plots
    make_power_plots(
      cell_types = unique(seurat_obj[[glue("RNA_snn_res.{res}")]])[[1]],
      outdir = pow_dir,
      markerfile_command = quote(glue("clust{cell_type}.csv")),
      marker_dir = d, plot_class = "avg_log2FC")

    # update list entry
    seurat_list[[i]] <- seurat_obj
  }
  return(seurat_list)
}

ds_list <- list("control" = embryoid_ctl, "nicotine" = embryoid_nic)
ds_list <- perform_seurat_analysis(ds_list, npcs, image_dir, obj_dir, csv_dir)

# save the outputs (override previous RDSs)
saveRDS(ds_list[1], file.path(data_dir, "embryoid_control2.rds"))
saveRDS(ds_list[2], file.path(data_dir, "embryoid_nicotine2.rds"))

# now subset FindMarkers values for p value less than 1%

saveRDS(ds_list, file.path(base_dir, "ds_list.rds"))

### integrate the two datasets

# # select features that are repeatedly variable across datasets for integration
# features <- SelectIntegrationFeatures(object.list = ds_list)
# anchors <- FindIntegrationAnchors(object.list = ds_list,
#                                   anchor.features = features)
# # create combined object
# embryoid <- IntegrateData(anchorset = anchors)

# # specify that we will perform downstream analysis on the corrected data note
# # that the original unmodified data still resides in the 'RNA' assay
# DefaultAssay(embryoid) <- "integrated"

# # Run processing on combined
# embryoid <- FindVariableFeatures(object = embryoid,
#                                  mean.function = ExpMean,
#                                  dispersion.function = LogVMR,
#                                  x.low.cutoff = min_avg_expr,
#                                  x.high.cutoff = max_avg_expr,
#                                  y.cutoff = min_dispersion)
# embryoid <- ScaleData(embryoid, features = rownames(embryoid),
#                       vars.to.regress = c("percent.mt", "nCount_RNA"))
# embryoid <- RunPCA(embryoid,
#                    features = VariableFeatures(embryoid),
#                    npcs = npcs)
# embryoid <- RunTSNE(embryoid, reduction = "pca", dims = 1:npcs)
# embryoid <- FindNeighbors(embryoid, reduction = "pca", dims = 1:npcs)
# embryoid <- FindClusters(embryoid, resolution = 0.8)

# markers <- FindMarkers(embryoid, logfc.threshold = 0.25, min.pct = 0.25,
#                        only.pos = TRUE, assay = "RNA")

# # Visualization
# p1 <- DimPlot(embryoid, reduction = "umap", group.by = "project")
# p2 <- DimPlot(emnbryooid, reduction = "umap", label = TRUE, repel = TRUE)
# p1 + p2

# # For performing differential expression after integration, we switch back to the original
# # data
# DefaultAssay(embryoid) <- "RNA"

# # create an outdir directory for FindMarkers run
# d <- file.path(csv_dir, "combined_markers_per_clust")

# # for cluster in whatever the combined cluster names are
# for (clust in embryoid$RNA_snn_res_0.8) {
#   clust_markers <- FindConservedMarkers(embryoid, ident.1 = clust,
#                                         grouping.var = "project",
#                                         verbose = FALSE)
#   head(clust_markers)

# }


# # label param set to cell type in featureplot
