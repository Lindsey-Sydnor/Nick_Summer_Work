#!/usr/bin/env Rscript

# Purpose:
#   TODO

# create_package.R first. Integration is personalized package.

# load packages
my_packages <- c("Seurat", "dyplr", "patchwork", "glue", "ggplot2", "clustree",
                 "ggraph", "pryr", "integration")

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

# setting default dirs
data_dir <- file.path(base_dir, "data", "organoid") # premade, stores raw data
# refs to output dirs
image_dir <- file.path(base_dir, "output", "images")
obj_dir <- file.path(base_dir, "output", "objs")
csv_dir <- file.path(base_dir, "output", "CSVs")
umap_dir <- file.path(image_dir, "umaps")

# creating output directories
dirs <- c(image_dir, obj_dir, csv_dir, umap_dir)
for (d in dirs){
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# presets
npcs <- 30

# load object -- Nick's / server's RDS was not working for me.
# organoid <- readRDS(file.path(data_dir, "organoid.rds"))

# command <- intersect(x = c("FindIntegrationAnchors", "FindTransferAnchors"),
#                            y = Command(object = organoid))[1]

# Command(object = organoid, command = command, value = "normalization.method")


# this Nick organoid RDS wasn't working for me... start with their RDS's and see
mg_rds <- readRDS(file.path(data_dir, "GSM4524697_NAY6153A1_125.rds"))
mg_rds <- UpdateSeuratObject(mg_rds)
mg_rds$stim <- "MG"
ctrl_rds <- readRDS(file.path(data_dir, "GSM4524699_NAY6153A2_678.rds"))
ctrl_rds <- UpdateSeuratObject(ctrl_rds)
ctrl_rds$stim <- "CTRL"

# for each, perform preprocessing
rdss <- list("MG_dataset" = mg_rds, "CTRL_dataset" = ctrl_rds)
for (i in seq_along(rdss)) {
  rdss[[i]] <- PercentageFeatureSet(rdss[[i]], pattern = "^MT-",
                                    col.name = "percent.mt")

  # make outdir
  d <- file.path(image_dir, "pre_integrated", "pre_QC", names(rdss[i]))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

  # create violin plots of features of interest per hash_ID
  to_test <- c("percent.mt", "nFeature_RNA", "nCount_RNA")
  for (test in to_test) {
    p <- VlnPlot(rdss[[i]], features = test)
    ggsave(plot = p, filename = file.path(d, glue("{test}.png")))
  }

  # subset same as in their analysis (no percent.mt considered):
  # remove non-singlet cells & threshold nFeature_RNA
  rdss[[i]] <- subset(rdss[[i]],
                      subset = nFeature_RNA > 500 &
                        rdss[[i]]$hto_classification_global == "Singlet")

  # retain active identity as sample-identifying antibody
  Idents(rdss[[i]]) <- "hash_ID"

  # make outdir
  d <- file.path(image_dir, "pre_integrated", "post_QC", names(rdss[i]))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

  # create violin plots of features of interest per hash_ID
  to_test <- c("percent.mt", "nFeature_RNA", "nCount_RNA")
  for (test in to_test) {
    p <- VlnPlot(rdss[[i]], features = test)
    ggsave(plot = p, filename = file.path(d, glue("{test}.png")))
  }

  # perform rest of preprocessing
  rdss[[i]] <- SCTransform(rdss[[i]], vars.to.regress = "percent.mt",
                           verbose = FALSE)
  rdss[[i]] <- RunPCA(rdss[[i]], verbose = FALSE)
  rdss[[i]] <- RunUMAP(rdss[[i]], dims = 1:npcs, verbose = FALSE)
  rdss[[i]] <- FindNeighbors(rdss[[i]], dims = 1:npcs, verbose = FALSE)
  rdss[[i]] <- FindClusters(rdss[[i]], resolution = 0.8, verbose = FALSE)
}

# unpack results
mg_rds <- rdss[[1]]
ctrl_rds <- rdss[[2]]

# create comparative plot
p1 <- DimPlot(mg_rds, label = TRUE, reduction = "umap") + NoLegend() +
  ggtitle("MG Dataset")
p2 <- DimPlot(ctrl_rds, label = TRUE, reduction = "umap") + NoLegend() +
  ggtitle("Unencapsulated Dataset")
p <- plot_grid(p1, p2)

# mkdir
d <- file.path(image_dir, "pre_integrated")
dir.create(d, showWarnings = FALSE)
# save
ggsave(filename = file.path(d, "comparative_umaps.png"), plot = p)

# prepare SCT assay for integration
features <- SelectIntegrationFeatures(object.list = rdss, nfeatures = 3000)
rdss <- PrepSCTIntegration(object.list = rdss, anchor.features = features)

# integrate datasets
rdss_anchors <- FindIntegrationAnchors(object.list = rdss,
                                       normalization.method = "SCT",
                                       anchor.features = features)
# create combined organoid object
organoid <- IntegrateData(anchorset = rdss_anchors,
                          normalization.method = "SCT")

# create integrated reductions
organoid <- RunPCA(organoid, verbose = FALSE)
organoid <- RunUMAP(organoid, reduction = "pca", dims = 1:npcs,
                    resolution = 0.8) # default res

# compare stim and clustering
p1 <- DimPlot(organoid, reduction = "umap", group.by = "stim")
p2 <- DimPlot(organoid, reduction = "umap", group.by = "SCT_snn_res.0.8",
              label = TRUE, repel = TRUE)
p <- p1 + p2

# create dir
d <- file.path(image_dir, "post_integrated", "misc")
dir.create(d, recursive = TRUE, showWarnings = FALSE)
ggsave(filename = file.path(d, "stim_v_res0.8.png"), plot = p, width = 12,
       height = 8, units = "in")

saveRDS(organoid, file = file.path(obj_dir, "correct_integrated_organoid.rds"))

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

# make outdir
misc_dir <- file.path(image_dir, "misc")
dir.create(misc_dir, showWarnings = FALSE)

# save elbowplot to misc dir
p <- ElbowPlot(organoid)
ggsave(filename = file.path(misc_dir, "elbow_plot.png"), plot = p)

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

# save checkpoint
saveRDS(object = organoid, file = file.path(obj_dir, "organoid_updated.rds"))
organoid <- readRDS(file.path(obj_dir, "organoid_updated.rds"))

# run clustree for ideal resolution info
clust_diagram <- integration::format_legend(clustree(organoid,
                  prefix = "integrated_snn_res."), text_size = 7)

# save to misc dir
ggsave(filename = file.path(misc_dir, "clustree.png"), plot = clust_diagram)

# rerun processing pipeline @ multiple resolutions on integrated assay:
resolutions <- seq(0.1, 1, by = 0.1)

# create outdir for FindMarkers runs
dir.create(file.path(csv_dir, "FindMarkers_run"), showWarnings = FALSE)

# create image outdir
d <- file.path(umap_dir, "pre_QC", "by_res")
dir.create(d, showWarnings = FALSE)

# make an info.txt specifying params for FindMarkers run
file_name <- "info.txt"
text <-
  "FindMarkers ran with these parameters:\
  \
    FindMarkers(organoid, ident.1 = clust, test.use = 'wilcox',\
                min.pct = 0.1, only.pos = FALSE, max.cells.per.ident = Inf,\
                verbose = FALSE)"

# write the text to the file
writeLines(text, file.path(csv_dir, "FindMarkers_run", file_name))

# Find clusters @ each resolution
for (r in resolutions) {
  # set keys and title per resolution
  cluster_key <- paste0("res", r)
  plot_title <- paste0("Resolution ", r)

  # run FindClusters @ each r
  organoid <- FindClusters(organoid, verbose = FALSE, resolution = r)
  organoid[[glue("integrated_snn_res.{r}")]] <- organoid$seurat_clusters

  # set active ident to r of newly-computed clusters
  Idents(organoid) <- glue("integrated_snn_res.{r}")

  # create umap @ each r
  p <- DimPlot(organoid, reduction = "umap",
               group.by = glue("integrated_snn_res.{r}")) + NoAxes()

  # save to "d" (set pre loop)
  ggsave(filename = file.path(d, glue("{r}.png")), plot = p)

  # make a directory for each res's FindMarkers output
  fm_dir <- file.path(csv_dir, "FindMarkers_run", glue("res{r}"))
  dir.create(fm_dir, showWarnings = FALSE)

  # set default assay to RNA for marker analysis
  DefaultAssay(organoid) <- "RNA"

  # run FindMarkers for each cluster
  for (clust in (levels(organoid))){
    print(glue("Running FindMarkers on cluster {clust}..."))

    # run FindMarkers on RNA assay
    cell_markers <- FindMarkers(organoid, assay = "RNA", ident.1 = clust,
                                test.use = "wilcox", min.pct = 0.1,
                                only.pos = FALSE,
                                max.cells.per.ident = Inf,
                                verbose = FALSE)
    # save to csv
    write.csv(x = cell_markers, file = file.path(d, glue("{clust}.csv")))
  }
}

DefaultAssay(VLMC_brain) <- "RNA" # should already be the case

for (res in resolutions) {
  print(res)
  cluster_key <- paste0("res", res)
  plot_title <- paste0("Resolution ", res)
  VLMC_brain <- FindNeighbors(VLMC_brain, dims = 1:npcs, verbose = FALSE)
  VLMC_brain <- FindClusters(VLMC_brain, verbose = FALSE, resolution = res) 
                              # name = cluster_key)
  res_col <- paste0("RNA_snn_res.", res)
  p <- DimPlot(object = VLMC_brain, reduction = "umap",
              group.by = res_col, pt.size = 0.5, label = TRUE)
  p <- p + theme(legend.position = "bottom") + NoAxes()
  # Save the modified ggplot object to file
  ggsave(filename = file.path(d, glue("{cluster_key}.png")), plot = p,
         width = 8, height = 8, dpi = 300)
}







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