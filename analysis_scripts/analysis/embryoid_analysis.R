#!/usr/bin/env Rscript

# Author: Lindsey Sydnor

# Purpose:
#   TODOs included...

# integration package is custom

# load packages
my_packages <- c("Seurat", "Matrix", "glue", "ggplot2", "integration",
                 "patchwork")

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

# presets
npcs <- 30

# set paths
base_dir <- getwd() # TODO: CHANGE FOR GIT
data_dir <- file.path(base_dir, "data", "embryoid") # premade, stores raw data
image_dir <- file.path(base_dir, "embryoid_output", "images")
csv_dir <- file.path(base_dir, "embryoid_output", "CSVs")
obj_dir <- file.path(base_dir, "embryoid_output", "objs")

# creating output directories
dirs <- c(image_dir, obj_dir, csv_dir)
for (d in dirs){
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# read in data
counts <- readMM(file.path(data_dir, "GSM3573649_D_matrix.mtx"))
# TODO: ^ large file, in .gitignore for now
barcodes <- read.table(file.path(data_dir, "GSM3573649_D_barcodes.tsv"))
genes <- read.table(file.path(data_dir, "GSM3573649_D_genes.tsv"))

# construct counts matrix for Seurat
colnames(counts) <- barcodes$V1
rownames(counts) <- genes$V2 # extract gene symbols, ignore ensembl IDs

# TODO: interesting... non-unique features
embryoid <- CreateSeuratObject(counts = counts, project = "embryoid",
                               min.cells = 3, min.features = 200)

embryoid <- PercentageFeatureSet(embryoid, pattern = "^MT-",
                                 col.name = "percent.mt")
pre_sub <- ncol(embryoid)
embryoid <- subset(embryoid, subset = percent.mt < 10)
post_sub <- ncol(embryoid)
print(glue("Loss of {pre_sub - post_sub} after percent.mt thresh of 10"))

# perform preprocessing
embryoid <- SCTransform(embryoid, vars.to.regress = "percent.mt",
                        verbose = FALSE)
embryoid <- RunPCA(embryoid, assay = "SCT", npcs = npcs)
embryoid <- RunUMAP(embryoid, dim = 1:npcs)
embryoid <- FindNeighbors(embryoid, dims = 1:npcs, assay = "SCT",
                          verbose = FALSE)

# germ layer gene markers
germ_markers <- list("Endoderm" = c("SOX17", "FOXA2", "CXCR4", "GATA4"),
                     "Mesoderm" = c("NCAM1", "TBXT"),
                     "Ectoderm" = c("NES", "PAX6"))

d <- file.path(image_dir, "umaps", "by_res")
dir.create(d, recursive = TRUE, showWarnings = FALSE)

# create outdir for germ-layer-specific analysis
germ_dir <- file.path(image_dir, "germ_layers")
dir.create(germ_dir, showWarnings = FALSE)

# run clustering and germ layer marker visualizations at multiple resolutions
resolutions <- seq(0.1, 1, by = 0.1)
for (r in resolutions) {
  # create outdir for germ-layer-specific analysis
  res_dir <- file.path(image_dir, "germ_layers", "expr_overlays",
                       glue("res{r}"))
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

  # set keys and title per resolution
  cluster_key <- paste0("res", r)

  # run FindClusters
  embryoid <- FindClusters(embryoid, verbose = FALSE, resolution = r)
  
  # set active ident to r of newly-computed clusters
  Idents(embryoid) <- glue("SCT_snn_res.{r}")

  # create UMAP
  p <- DimPlot(embryoid, reduction = "umap",
                group.by = glue("SCT_snn_res.{r}"), label = TRUE) + NoAxes() +
    NoLegend() +
    ggtitle(glue("Resolution {r}"))
  
  ggsave(filename = file.path(d, glue("{r}.png")), plot = p1)

  p1 <- DimPlot(embryoid, reduction = "umap",
              group.by = glue("SCT_snn_res.{r}"), label = TRUE) + NoAxes() +
    NoLegend()

  # save to "d" (set pre loop)
  ggsave(filename = file.path(d, glue("{r}.png")), plot = p1)

  # germ marker exploration
  for (germ_layer in names(germ_markers)) {
    # set shared title for magick_overlay plots
    shared_title <- glue("Resolution {r} {cell_type} Canonical Gene",
                         " Expression")

    # make an overlay with a FeaturePlot showing germ layer markers
    p <- FeaturePlot(embryoid, features = germ_markers[[germ_layer]]) +
      plot_annotation(glue("{names(germ_markers[germ_layer])} Markers at",
                           " {shared_title}"),
                      theme = theme(plot.title = element_text(hjust = 0.5)))

    # add extra width if two plots in patchwork for better visual
    if (p$patches$layout$ncol == 2) {
      ggsave(filename = file.path(germ_dir, glue("{germ_layer}_feature.png")),
             plot = p, width = 10, height = 6, units = "in")
    } else {
      # save FeaturePlot as is
      ggsave(filename = file.path(germ_dir, glue("{germ_layer}_feature.png")),
             plot = p)
    }

    p1 <- p1 + ggtitle(shared_title)
    
    # get ranges of plot for expression overlays
    xrange <- layer_scales(p1)$x$range$range
    yrange <- layer_scales(p1)$y$range$range

    ### make UMAP of all cells expressing this germ layer's canonical genes
    # determine which genes are present in embryoids
    matching_genes <- intersect(rownames(embryoid), germ_markers[[germ_layer]])
    # include any cells expressing any canonical genes above threshold:
    threshold <- 0
    # initialize
    selected_cells <- c()
    # for every germ layer canonical gene found in the dataset...
    for (matched in matching_genes) {
      before <- length(selected_cells)
      selected_cells <- c(selected_cells,
                colnames(embryoid)[embryoid@assays$RNA@data[matched, ] >
                         threshold])
      selected_cells <- unique(selected_cells) # prune cells already on list
      after <- length(selected_cells)
      print(glue("Num cells containing gene {matched}: {after - before}"))
    }
    if (length(selected_cells) > 0) {
      # Create a subset containing only the selected cells
      subset_obj <- embryoid[, selected_cells]
      # plot
      p2 <- DimPlot(subset_obj, reduction = "umap", label = FALSE) +
        NoLegend() + ggtitle(shared_title) +
        xlim(xrange) + ylim(yrange) +
        theme(plot.title = element_text(hjust = 0.5)) +
        NoAxes()
      im <- integration::magick_overlay(main_plt = p2, trans_plt = p1,
                                          dest = res_dir,
                                          filename = glue("{germ_layer}.png"),
                                          x_dim = 14, y_dim = 14)
    }
  }

  # # make a directory for each res's FindMarkers output
  # fm_dir <- file.path(csv_dir, "FindMarkers_run", cluster_key)
  # dir.create(fm_dir, recursive = TRUE, showWarnings = FALSE)

  # for (clust in levels(embryoid)) {
  #   markers <- FindMarkers(embryoid, assay = "SCT", ident.1 = clust,
  #                          verbose = FALSE)
  #   write.csv(x = markers, file = file.path(fm_dir, glue("{clust}.csv")))
  # }
}

# save checkpoint for iPSC_germ_markers.R
saveRDS(embryoid, file.path(obj_dir, "embryoid_update.rds"))

# has more cells excluded from analysis than their thresholding standards
embryoid <- readRDS(file.path(obj_dir, "embryoid_update.rds"))


### Cell-type-specific analysis
# canonical cell type markers extracted from their paper
cell_markers <- list("Liver_Progenitors" = c("FRZB", "PTN"),
                     "Neural_Progenitors" = c("LHX2", "NR2F1"),
                     "Stromal_Progenitors" = c("SFRP2", "COL2A1"),
                     "Endothelial_Progenitors" = c("GDF15", "DDIT3"),
                     "Epithelial_Progenitors" = c("IGFBP5", "PODXL"),
                     "Muscle_Progenitors" = c("HAPLN1", "S100A11"),
                     "Undifferentiated_Stem-Like_Cells" = c("TERF1", "POU5F1"),
                     "Undetermined_Cells" = c("ACTB", "TUBB"))

submarkers <- list("Forebrain_Development" = c("LHX5", "HESX1"),
                   "Proliferating_Neural_Stem_Cells" = c("HMGB2", "PTTG1"),
                   "Sensory_Neuron_Development" = c("HNRNPH1", "PTPRS"),
                   "Eye_Development" = c("LHX2", "NR2F1"),
                   "Smooth_Muscle" = c("HAPLN1", "S100A11"),
                   "Cardiac_Muscle" = c("ZFHX3", "NR2F2"),
                   "Eye_Epithelial" = c("IGFBP5", "HES1"),
                   "Stem-Like_Epithelial" = c("B3GNT7", "PODXL"))

# create outdir for cell-type-specific analysis
cell_dir <- file.path(image_dir, "cell_markers")
dir.create(cell_dir, showWarnings = FALSE)

# create outdir for cell-type-specific analysis
sub_dir <- file.path(image_dir, "submarkers")
dir.create(sub_dir, showWarnings = FALSE)

for (cell_type in names(cell_markers)) {
    # make an overlay with a FeaturePlot showing germ layer markers
    p <- FeaturePlot(embryoid, features = cell_markers[[cell_type]]) +
      plot_annotation(glue("{gsub('_', ' ', names(cell_markers[cell_type]))}",
                           " Canonical Markers"),
                      theme = theme(plot.title = element_text(hjust = 0.5)))

    # add extra width if two plots in patchwork for better visual
    if (p$patches$layout$ncol == 2) {
      ggsave(filename = file.path(cell_dir,
             glue("{cell_type}_feature.png")),
             plot = p, width = 10, height = 6, units = "in")
    } else {
      # save FeaturePlot as is
      ggsave(filename = file.path(cell_dir,
             glue("{cell_type}_feature.png")),
             plot = p)
    }
}

for (cell_type in names(submarkers)) {
    # make an overlay with a FeaturePlot showing germ layer markers
    p <- FeaturePlot(embryoid, features = submarkers[[cell_type]]) +
      plot_annotation(glue("{gsub('_', ' ', names(submarkers[cell_type]))} ",
                           "Canonical Markers"),
                      theme = theme(plot.title = element_text(hjust = 0.5)))

    # add extra width if two plots in patchwork for better visual
    if (p$patches$layout$ncol == 2) {
      ggsave(filename = file.path(sub_dir, glue("{cell_type}_feature.png")),
             plot = p, width = 10, height = 6, units = "in")
    } else {
      # save FeaturePlot as is
      ggsave(filename = file.path(sub_dir, glue("{cell_type}_feature.png")),
             plot = p)
    }
}

# run clustering and germ layer marker visualizations at multiple resolutions
resolutions <- seq(0.1, 1, by = 0.1)
for (r in resolutions) {
  # create outdir for germ-layer-specific analysis
  res_dir <- file.path(image_dir, "cell_types", "expr_overlays",
                       glue("res{r}"))
  dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

  # set keys and title per resolution
  cluster_key <- paste0("res", r)

  # run FindClusters
  # embryoid <- FindClusters(embryoid, verbose = FALSE, resolution = r)
  
  # set active ident to r of newly-computed clusters
  Idents(embryoid) <- glue("SCT_snn_res.{r}")

  # create UMAP

  p1 <- DimPlot(embryoid, reduction = "umap",
              group.by = glue("SCT_snn_res.{r}"), label = TRUE) + NoAxes() +
    NoLegend()

  # germ marker exploration
  for (cell_type in names(cell_markers)) {
    # set shared title for magick_overlay plots
    shared_title <- glue("Resolution {r} {cell_type} Canonical Gene",
                         " Expression")

    p1 <- p1 + ggtitle(shared_title)
    
    # get ranges of plot for expression overlays
    xrange <- layer_scales(p1)$x$range$range
    yrange <- layer_scales(p1)$y$range$range

    ### make UMAP of all cells expressing this germ layer's canonical genes
    # determine which genes are present in embryoids
    matching_genes <- intersect(rownames(embryoid), cell_markers[[cell_type]])
    # include any cells expressing any canonical genes above threshold:
    threshold <- 0
    # initialize
    selected_cells <- c()
    # for every germ layer canonical gene found in the dataset...
    for (matched in matching_genes) {
      before <- length(selected_cells)
      selected_cells <- c(selected_cells,
                      colnames(embryoid)[embryoid@assays$RNA@data[matched, ] >
                                          threshold])
      selected_cells <- unique(selected_cells) # prune cells already on list
      after <- length(selected_cells)
      print(glue("Num cells containing gene {matched}: {after - before}"))
    }
    if (length(selected_cells) > 0) {
      # Create a subset containing only the selected cells
      subset_obj <- embryoid[, selected_cells]
      # plot
      p2 <- DimPlot(subset_obj, reduction = "umap", label = FALSE) +
        NoLegend() + ggtitle(shared_title) +
        xlim(xrange) + ylim(yrange) +
        theme(plot.title = element_text(hjust = 0.5)) +
        NoAxes()
      im <- integration::magick_overlay(main_plt = p2, trans_plt = p1,
                                        dest = res_dir,
                                        filename = glue("{cell_type}.png"),
                                        x_dim = 14, y_dim = 14)
    }
  }

  # # make a directory for each res's FindMarkers output
  # fm_dir <- file.path(csv_dir, "FindMarkers_run", cluster_key)
  # dir.create(fm_dir, recursive = TRUE, showWarnings = FALSE)

  # for (clust in levels(embryoid)) {
  #   markers <- FindMarkers(embryoid, assay = "SCT", ident.1 = clust,
  #                          verbose = FALSE)
  #   write.csv(x = markers, file = file.path(fm_dir, glue("{clust}.csv")))
  # }
}

# # IDK WHERE TO PUT VIOLIN PLOTS YET
# # Visualize canonical marker genes as violin plots.
# p <- VlnPlot(embryoid, features = germ_markers$Endoderm, pt.size = 0.2,
#               ncol = 2) + ggtitle("Resolution {r} Endoderm Markers")

# p <- VlnPlot(embryoid, features = germ_markers$Mesoderm, pt.size = 0.2,
#               ncol = 2) + ggtitle("Resolution {r} Mesoderm Markers")

# p <- VlnPlot(embryoid, features = germ_markers$Ectoderm, pt.size = 0.2,
#               ncol = 2) + ggtitle("Resolution {r} Ectoderm Markers")

# d <- file.path(image_dir, "germ_markers", "violin_plots", glue("res{r}"))
# dir.create(d, showWarnings = FALSE)
# # ggsave(file = file.path(d, "endoderm"),


# # Visualize canonical marker genes as violin plots.
# p <- VlnPlot(embryoid, features = germ_markers$Endoderm, pt.size = 0.2,
#              ncol = 2)

# d <- file.path(image_dir, "germ_markers", "violin_plots", glue("res{r}"))
# dir.create(d, showWarnings = FALSE)