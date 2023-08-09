# Contains:
#   - save_fxn
#   - format_legend
#   - magick_overlay
#   - query_integrated_reference
#   - runBIC
#   - run_qc
#   - findmarkers_gene2cell_mapping
#   - make_expression_overlays
# TODO: 
#     Flush out descriptions. Format file so linter stops yelling.

#' save_fxn
#'
#' @description Save a ggplot object to the appropriate location(s)
#' @param obj (ggplot obj): \cr
#'  ggplot object to save
#' @param file (string): \cr
#'  file name to save ggplot as (include file extension: "png" only)
#' @param w (int or double): \cr
#'  width of PNG to save in inches
#' @param h (int or double): \cr
#'  height of PNG to save in inches
#' @param dirs (string or character vector): \cr
#'  directory (or directories) to save the plot in
#' @import ggplot2 glue
#' @return N/A
#' @examples
#'  save_fxn(main_plt, file = "tmp_plt1.png", dirs = c("~/endothelial/images",
#'         "~/examples/images"))
#' @export
save_fxn <- function(obj, file, dirs, ..., w = 8, h = 5) {
  for (d in dirs) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
    ggsave(plot = obj, filename = file.path(path = d, file), device = "png",
           units = "in", width = w, height = h)
  }
  # save 1, copy over rest (faster)
  if (length(dirs) > 1) {
    for (i in 2:length(dirs)) {
      file.copy(from = paste0(dirs[1], file),
                to = paste0(dirs[i], file), overwrite = TRUE)
    }
  }
}

#' format_legend
#'
#' @description Ensure single-column, appropriately-sized legend in ggplot
#' @param my_plot (ggplot obj): \cr
#'  ggplot object to return
#' @param point_size (float): \cr
#'  size to make points in legend
#' @param text_size (int or double): \cr
#'  size to make text in legend
#' @param space_legend (float): \cr
#'  linespacing of legend
#' @import ggplot2
#' @return ggplot object
#' @examples 
#' p1 <- format_legend(DimPlot(seurat_obj, reduction = "umap",
#'                    group.by = "sample", label = TRUE, label.size = 2,
#'                    repel = TRUE))
#' @export
format_legend <- function(my_plot, ..., point_size = 0.5, text_size = 5,
                          space_legend = point_size * 0.5) {
  my_plot <- my_plot +
    guides(shape = guide_legend(override.aes = list(size = point_size)),
           color = guide_legend(override.aes = list(size = point_size),
                                ncol = 1), # also makes legend a single column
           fill = guide_legend(override.aes =  list(ncol = 1))) +
    theme(legend.title = element_text(size = text_size),
          legend.text  = element_text(size = text_size),
          legend.key.size = unit(space_legend, "lines"))
  return(my_plot)
}

#' magick_overlay
#'
#' @description Create an overlay of two plots, with `transparency`%
#'  translucency, to visualize both simultaneously.
#' \cr
#' NOTE: Ideally, make ggtitle the same value on the two plots to overlay or
#'       nonexistent. Additionally, best to supply plots without axes unless
#'       identical prior to function call.
#' @param trans_plt (ggplot obj): \cr ggplot to display with `transparency`%
#'  translucency
#' @param main_plt (ggplot obj): \cr ggplot to display with 0% translucency
#' @param dest (string or character vector): \cr directory or directories to 
#'  save overlayed magick image object (will create directories if not yet in
#'  existence)
#' @param filename (string): \cr file name to save overlay plot to (include file
#'  extension: "png" only)
#' @param x_dim (int or double = 14): \cr positive x axis limit for overlay
#'  plot. \cr
#'  Example: \cr
#'    To create a plot with x values from \[-17, 17\], supply: `x_dim` = 17.
#'  Assumed symmetric about the origin.
#' @param y_dim (int or double = 17): \cr positive y axis limit for overlay
#'  plot. \cr
#'  Example: \cr 
#'    To create a plot with y values from \[-17, 17\], supply `y_dim` = 17.
#'  Assumed symmetric about the origin.
#' @param transparency (int or double = 20): \cr percent transparency to make
#'  the transparent plot in the overlay
#' @param return_im (bool = FALSE) whether to return magick image object
#' @import ggplot2 glue magick
#' @return magick image object of overlayed plots (`main_plt` + `trans_plt`) if
#'  `return_im` == TRUE. Else N/A
#' @examples
#' im <- magick_overlay(main_plt = p1, trans_plt = p2,
#'                      dest = "~/endothelial/images/overlays/by_age",
#'                      filename = "10_PCW.png")
#' @export
magick_overlay <- function(main_plt, trans_plt, dest, filename, ..., x_dim = 14,
                           y_dim = 17, transparency = 20, return_im = FALSE) {
  # assummed symmetric about x and y axes
  # impose same dimensions, even if they were already included
  suppressMessages(main_plt + xlim(-x_dim, x_dim) + ylim(-y_dim, y_dim))
  suppressMessages(trans_plt + xlim(-x_dim, x_dim) + ylim(-y_dim, y_dim))
  # save then reload using magick, as no fxn for ggplot -> magick-image type
  save_fxn(main_plt, dirs = dest, file = "tmp_plt1.png")
  save_fxn(trans_plt, dirs = dest, file = "tmp_plt2.png")
  main_plt <- image_read(paste0(dest, "/tmp_plt1.png"))
  trans_plt <- image_read(paste0(dest, "/tmp_plt2.png"))
  img <- image_join(main_plt, trans_plt)
  img <- image_scale(img, "100X100")
  im <- image_composite(main_plt, composite_image = trans_plt,
                        operator = "dissolve", compose_args = 
                          glue("{toString(transparency)}%"),
                        gravity = "center")
  for (d in dest) {
    dir.create(d, showWarnings = FALSE, recursive = TRUE)
  }
  image_write(im, path = paste0(dest[1], "/", filename), format = "png")
  if (length(dest) > 1) {
    # copy the rest (faster than remake)
    for (i in 2:length(dest)) {
      file.copy(from = paste0(dest[1], "/", filename),
                to = paste0(dest[i], "/", filename), overwrite = TRUE)
    }
  }
  if (return_im == TRUE) {
    return(im)
  }
}

#' query_integrated_reference
#'
#' @description Integrate datasets (including or excluding `query` dataset) to
#'  create a reference against which to assign cell type labels to the `query`
#'  dataset. Datasets must be Seurat objects. Cell type labels must be under
#'  slot "cell_type" for all datasets. \cr
#'  Creates three directories (two subdirs of the first):  \cr
#'    - `out_dir`  \cr
#'    - `out_dir`/images  \cr
#'      dir to hold plots created by this process. These plots are designed
#'      to assess successful completion of the integration. \cr
#'    - `out_dir`/integration_objs
#'      dir to hold all integration objects created in the process, in case
#'      any step must be repeated or object(s) exported.
#' NOTES: Should have already run pca and other items on each object.
#' @param references (Seurat obj or list of Seurat objs): \cr dataset(s) to
#'  integrate. \cr
#'  NOTE: ensure `references` are normalized and variable features are found.
#'  Only run \link[Seurat]{FindVariableFeatures} on unnormalized data.
#' @param query (Seurat obj): \cr  query dataset of cells to label
#' @param out_dir (str): \cr parent directory name to hold images and
#'  integration objects created
#' @param npcs (int or double = 30): \cr number of principal components with
#'  which to run \link[Seurat]{FindIntegrationAnchors},
#'  \link[Seurat]{IntegrateData}, \link[Seurat]{RunUMAP},
#'  \link[Seurat]{FindTransferAnchors}, and \link[Seurat]{TransferData}.
#' @param base_dir (string = \link{getwd()} (currrent working dir)
#'  or char vect): \cr directory or directories) to save the plot in
#' @param include_query (bool = TRUE): \cr whether or not to include the `query`
#'  dataset when building the reference atlas of cells and cell type labels to
#'  compare `query` against. \cr
#'  NOTE: Excluding `query`` will not allow your newly-assigned cell type
#'  labels to include what they were originally, even if they would be found as
#'  the ideal match by this workflow. All original cell type labels will be
#'  replaced.
#' @param assay (NULL or str): \cr find integration anchors for specified Seurat
#'  object assay. If not provided, will use DefaultAssay of each reference (&
#'  query, if include_query = TRUE) to integrate along.
#' @param return (bool = TRUE): \cr if return = TRUE, return the query and
#'  integrated Seurat objects created by this analysis.
#' @import Seurat
#' @return (named list): \cr list of named Seurat objects: \cr 
#' query =  query dataset with newly-constructed cell type labels in assay
#' predicted.cell_type) \cr
#' integrated dataset
#' @examples
#' output <- query_integrated_reference(references = list(seurat_obj1,
#'            seurat_obj2), query = seurat_obj3, npcs = 50, base_dir = "~/",
#'            include_query = FALSE, assay = "RNA", return = TRUE)
#' query <- output$query
#' integrated <- output$integrated
#' @export
query_integrated_reference <- function(references, query, out_dir, npcs = 30,
                                       base_dir = getwd(), include_query = TRUE,
                                       return = TRUE, assay = NULL,
                                       reference_assay = "RNA",
                                       query_assay = "RNA", scale = FALSE,
                                       num_anchor_feats = 2000) {
  # ADD DOCS FOR THESE:
  # assay e.g. c("RNA", "RNA") (assay for int one per object in object_list)
  # scale = FALSE (assumed ScaleData already ran)
  # num_anchor_feats (num of features to use as anchors for integrating
  #                   datasets)
  # 
  # if single reference, encapsulate it as list class for iteration
  if (class(references) == "Seurat") {
    references <- c(references)
  }
  for (ref in references) {
    DefaultAssay(ref) <- reference_assay
  }
  DefaultAssay(query) <- query_assay
  my_dims <- as.numeric(1:npcs)
  # create dirs as needed
  for (loc in c("images", "integration_objs", "tables")) {
    dir.create(file.path(base_dir, out_dir, loc), showWarnings = FALSE,
               recursive = TRUE)
  }
  # Table of orig cell_type proportions
  table1 <- table(query@meta.data$cell_type)
  write.csv(table1, file = file.path(base_dir, out_dir, "tables",
                                     "original_cell_type_numbers.csv"))
  # First map old UMAP cell_type labels
  p1 <- DimPlot(query, reduction = "umap", group.by = "cell_type",
                label = TRUE, repel = TRUE) + NoAxes()
  save_fxn(obj = p1, dirs = file.path(base_dir, out_dir, "images"),
           file = "old_labels_umap.png")
  message("Finding Integration Anchors...")
  if (include_query == TRUE) {
    # integration_anchors <- FindIntegrationAnchors(object.list =
    #                                             c(references, query),
    #                                             dims = my_dims)
    integration_anchors <- FindIntegrationAnchors(object.list =
                                                    c(references, query),
                                                  dims = my_dims, assay = assay,
                                                  verbose = TRUE, scale = scale,
                                                  anchor.features = num_anchor_feats)
  } else {
    integration_anchors <- FindIntegrationAnchors(object.list = references,
                                                  dims = my_dims)
    # rethink for single reference use case              
  }
  message("Saving Integration Anchors...")
  saveRDS(integration_anchors, file.path(base_dir, out_dir,
                                         "integration_objs", "integration_anchors.rds"))
  message("Integrating Datasets...")
  integrated <- IntegrateData(anchorset = integration_anchors, dims = my_dims,
                              verbose = TRUE)
  DefaultAssay(integrated) <- "integrated"
  # Set default to integrated assay. Variable features of this assay are
  # automatically determined / set during IntegrateData
  message("Running Standard Processing...")
  integrated <- ScaleData(integrated, verbose = FALSE)
  integrated <- RunPCA(integrated, npcs = npcs, verbose = FALSE)
  integrated <- RunUMAP(integrated, dims = my_dims, reduction = "pca", 
                        verbose = FALSE, return.model = TRUE)
  message("Saving Integrated Object...")
  saveRDS(integrated, file.path(base_dir, out_dir, "integration_objs",
                                "integrated.rds"))
  message("Plotting Integrated UMAP")
  p1 <- DimPlot(integrated, reduction = "umap", group.by = "cell_type",
                label = TRUE, repel = TRUE) + NoAxes()
  save_fxn(obj = p1, dirs = file.path(base_dir, out_dir, "images"),
           file = "integrated_umap.png")
  message("Finding Transfer Anchors...")
  # Transfer old cell type labels from reference dataset onto new query dataset
  transfer_anchors <- FindTransferAnchors(reference = integrated, query = query,
                                          dims = my_dims, reference.reduction = "pca", 
                                          reference.assay = reference_assay,
                                          query.assay = query_assay)
  message("Saving Transfer Anchors...")
  saveRDS(transfer_anchors, file.path(base_dir, out_dir, "integration_objs",
                                      "transfer_anchors.rds"))
  # slot(transfer_anchors, name = 'command')$reference.assay # integrated
  # slotNames(integration_anchors)
  # head(integration_anchors@anchors, n=2)
  # > length(integration_anchors@anchor.features)
  # [1] 1803
  # > dim(integration_anchors@anchors) # should equal num cells in smaller dataset
  # [1] 6998    5
  message("Mapping to Query...")
  query <- MapQuery(anchorset = transfer_anchors, reference = integrated,
                    query = query, refdata = list(cell_type = "cell_type"),
                    reference.reduction = "pca", reduction.model = "umap")
  saveRDS(query, file.path(base_dir, out_dir, "integration_objs",
                           "updated_query.rds"))       
  message("Plotting Query Projected on Integrated Reference Atlas")
  p1 <- DimPlot(integrated, reduction = "umap", group.by = "cell_type", 
                label = TRUE, label.size = 3, repel = TRUE) + NoLegend() + 
    NoAxes() + ggtitle("Reference Annotations")
  p2 <- DimPlot(query, reduction = "ref.umap", group.by = "predicted.cell_type",
                label = TRUE, label.size = 3, repel = TRUE) + NoLegend() +
    NoAxes() + ggtitle("Query Transferred Labels")
  p3 <- p1 + p2
  save_fxn(p3, dirs = file.path(base_dir, out_dir, "images"),
           file = "umaps_compared.png")
  # map predicted.cell_types
  p1 <- DimPlot(query, reduction = "umap", group.by = "predicted.cell_type",
                label = TRUE, repel = TRUE) + NoAxes()
  save_fxn(obj = p1, dirs = file.path(base_dir, out_dir, "images"),
           file = "new_labels_umap.png")
  # table of predicted cell_type proportions
  table1 <- table(query@meta.data$predicted.cell_type)
  write.csv(table1, file = file.path(base_dir, out_dir, "tables",
                                     "predicted_cell_type_numbers.csv"))
  # write out previous and new cell_type numbers
  if (return == TRUE) {
    return(list(query = query, integrated = integrated))
  }
}

#' runBIC
#' @description
#'  Find the ideal number of variable features to use for a Seurat object's
#'  downstream analysis. \cr
#'  Process: Run \link[Seurat]{FindVariableFeatures} to find `nfeature_max`
#'  number of variable features. Then for each number of variable features: 
#'  \link[Seurat]{ScaleData}, run \link[Seurat]{SCTransform} on the variable
#'  features, run PCA, find the proportion of explained variance of the total
#'  dataset's variance for the first `npcs` components, calculate total variance
#'  of the dataset's "scale.data" slot, and then calculate the BIC value for
#'  that total PCA model. Compare BIC values across each number of variable
#'  features. The number of variable features corresponding to the lowest BIC
#'  value resulted in the "best" BIC model. If a plateau of BIC values is
#'  reached, the lowest number of variable features will be preferred. A plot
#'  of bic value per number of variable features will be produced in
#'  `file_path`.
#' @param obj (Seurat object): \cr
#'  Seurat object on which to run BIC analysis.
#' @param nfeature_min (int or double = 1500): \cr
#'  Minimum number of variable features to test.
#' @param nfeature_max (int or double = 5000): \cr
#'  Maximum number of variable features to test.
#' @param increment (int or double = 500): \cr
#'  Test ideal VariableFeature nfeature number between nfeature_min and
#'  nfeature_max by increments of `increment`
#' @param regress (NULL, string, or char vector): \cr
#'  Vector of variables to regress in \link[Seurat]{ScaleData} call, if any.
#' @param file_path (string = "~/BIC.png"): \cr
#'  Path to save BIC comparison png
#' @import Seurat glue ggplot2
#' @return (named list) (optimal_nfeatures, bic) \cr
#'  optimal_nfeatures: Ideal number of variable features as calculated via BIC.
#'  \cr bic: vector of bic values for each number of variable features tested.
#' @examples
#'  nvar_features <- runBIC(seurat_obj, file_path = glue("/active/aldinger_k/",
#'                          "kimslab_temp/scRNA-seq/brain-vasc/vascular-dev/",
#'                          "lsyd/test_save.png")
#' @export
runBIC <- function(obj, ..., nfeature_min = 1500, nfeature_max = 5000,
                   increment = 500, regress = NULL, plt = TRUE, npcs = 50,
                   file_path = "~/BIC.png") {
  DefaultAssay(obj) <- "RNA" # SCTransform won't run on integrated assay
  obj <- FindVariableFeatures(obj, nfeatures = nfeature_max)
  obj <- ScaleData(obj, features = rownames(obj), vars.to.regress = regress)
  nvals <- (nfeature_max - nfeature_min) / increment + 1
  # Initialize
  bic <- integer(nvals)
  for (i in 1:nvals) {
    nfeatures <- nfeature_min + increment * (i - 1)
    message(glue("Running {nfeatures} features..."))
    scale_data_df <- as.data.frame(obj@assays$RNA@scale.data)
    rownames_sorted <- rownames(scale_data_df)[order(apply(scale_data_df, 1,
                                                           max), decreasing = TRUE)][1:nfeatures]
    # Get the row indices corresponding to the sorted row names
    message("Finding indices...")
    row_indices <- match(rownames_sorted, rownames(scale_data_df))
    # Use the row indices to extract the variable features
    var_features <- obj@assays$RNA@var.features[row_indices[1:nfeatures]]
    message("Running SCTransform...")
    seurat_subset <- SCTransform(obj, residual.features = var_features,
                                 verbose = FALSE)
    # Take just pca reduction
    message("Running PCA...")
    pca <- RunPCA(seurat_subset, verbose = FALSE, npcs = npcs)[["pca"]]
    num_genes <- nrow(seurat_subset@assays$RNA@scale.data)
    num_cells <- ncol(seurat_subset@assays$RNA@scale.data)
    # Calc total variance of the data
    var_total <- sum(var(scale_data_df))
    num_components <- ncol(pca)
    # Calc the BIC value using the explained variance and num of components
    eig_values <- (pca@stdev)^2  # use eigen vals
    var_explained <- eig_values / var_total
    bic_value <- num_genes * num_cells * log(1 -
                                               sum(var_explained[1:num_components])) +
      num_components * log(num_genes * num_cells)
    # Store the BIC value
    bic[i] <- bic_value
  }
  # Plotting
  nfeature_vec = seq(nfeature_min, nfeature_max, by = increment)
  optimal_nfeatures <- nfeature_min + increment * (which.min(bic) - 1)
  if (plt) {
    message("Plotting")
    # Find the optimal number of features by the minimum BIC value
    
    # Do I want to output the optimal features dataset? Maybe make a flag
    # optimal_features <- seurat_obj@assays$RNA@var.features[1:optimal_nfeatures]
    x <- nfeature_vec
    y <- bic
    # p1 <- ggplot() +
    #   geom_point(aes(x = x, y = y), na.rm = TRUE) + 
    #   xlim(nfeature_min, nfeature_max) +
    #   scale_y_continuous(limits = c(min(y), max(y)))
    xlim <- c(nfeature_min, nfeature_max)
    ylim <- c(min(bic), max(bic))
    p1 <- ggplot() +
      geom_point(aes(x = nfeature_vec, y = bic)) +
      geom_vline(xintercept = optimal_nfeatures, linetype = "dashed",
                 color = "red") +
      scale_x_continuous(name = "Number of Variable Features", limits = xlim) +
      scale_y_continuous(name = "BIC Value", limits = ylim)
    ggsave(filename = file_path, plot = p1, 
           width = 6, height = 4, device = "png", dpi = 300)
  }
  return(list(optimal_nfeatures = optimal_nfeatures, bic = bic))
}

#' run_qc
#' @import Seurat reticulate ggplot2 tidyverse
#' @export
# adapted from online resource
# vln_by (what your groupings should be (i.e., by cluster of specific res,
# whatever))
run_qc <- function(seurat_obj, base_dir = getwd(), vln_by,
                   out_dir = "qc_outdir", mito_meta_name = "percent_mito") {
  dir.create(file.path(base_dir, out_dir), recursive = TRUE,
             showWarnings = FALSE)
  # Get Cell Cycle scoring
  s.genes <- cc.genes.updated.2019$s.genes # markers of s phase
  g2m.genes <- cc.genes.updated.2019$g2m.genes # markers of g2 and m phase
  seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes,
                                 g2m.features = g2m.genes,
                                 set.ident = FALSE)
  
  # create tibble for analysis
  qc.metrics <- as_tibble(seurat_obj[[]], rownames = "Cell.Barcode")
  p <- qc.metrics %>%
    arrange(seurat_obj[[mito_meta_name]]) %>%
    ggplot(aes(nCount_RNA, nFeature_RNA, colour = get(mito_meta_name))) +
    geom_point() +
    scale_color_gradientn(colors = c("black", "blue", "green2", "red",
                                     "yellow")) +
    ggtitle("QC metrics") +
    geom_hline(yintercept = 750) +
    geom_hline(yintercept = 2000) +
    labs(color = mito_meta_name) # set legend title
  
  # create output directory
  d <- file.path(base_dir, out_dir, "qc_metrics")
  dir.create(d, showWarnings = FALSE)
  
  ggsave(p, filename = file.path(base_dir, out_dir, "qc_metrics",
                                 "QC.png"), device = "png")
  # replot on log scale
  p <- qc.metrics %>%
    arrange(mito_meta_name) %>%
    ggplot(aes(nCount_RNA, nFeature_RNA, colour = get(mito_meta_name))) +
    geom_point(size = 0.7) +
    scale_color_gradientn(colors=c("black", "blue", "green2", "red", "yellow")) +
    ggtitle("Log-Scale QC metrics") +
    geom_hline(yintercept = 750) +
    geom_hline(yintercept = 2000) +
    scale_x_log10() + scale_y_log10() +
    labs(color = mito_meta_name)
  
  ggsave(p, filename = file.path(base_dir, out_dir, "qc_metrics",
                                 "log_QC.png"), device = "png")
  
  # plot complexity (UMIs / cell)
  qc.metrics <- qc.metrics %>% mutate(complexity = log10(nFeature_RNA) /
                                        log10(nCount_RNA))
  
  complexity.lm <- lm(log10(qc.metrics$nFeature_RNA)~log10(qc.metrics$nCount_RNA))
  qc.metrics <- qc.metrics %>% mutate(complexity_diff = log10(nFeature_RNA) -
                                        ((log10(qc.metrics$nCount_RNA) *
                                            complexity.lm$coefficients[2]) +
                                           complexity.lm$coefficients[1]))
  p <- qc.metrics %>% ggplot(aes(x = complexity_diff)) +
    geom_density(fill = "yellow")
  ggsave(p, filename = file.path(base_dir, out_dir, "qc_metrics",    
                                 "complex_v_dens.png"), device = "png")
  
  # add complexity to log qc mapping
  complexity_scale <- min(c(max(qc.metrics$complexity_diff),
                            0 - min(qc.metrics$complexity_diff)))
  p <- qc.metrics %>%
    mutate(complexity_diff = replace(complexity_diff,
                                     complexity_diff< -0.1,-0.1)) %>%
    ggplot(aes(x = log10(nCount_RNA), y = log10(nFeature_RNA),
               colour = complexity_diff)) +
    geom_point(size = 0.5) +
    geom_abline(slope = complexity.lm$coefficients[2],
                intercept = complexity.lm$coefficients[1]) +
    scale_colour_gradient2(low = "blue2", mid = "grey", high = "red2")
  ggsave(p, filename = file.path(base_dir, out_dir, "qc_metrics",    
                                 "log_QC_w_dens.png"), device = "png")
  
  # Cell Cyle Info
  p <- qc.metrics %>%
    ggplot(aes(x = S.Score, y = G2M.Score, color = Phase)) +
    geom_point() + ggtitle("Cell Cycle")
  ggsave(p, filename = file.path(base_dir, out_dir, "qc_metrics",
                                 "cell_cycle.png"), device = "png")
  
  # VLN plots
  Idents(seurat_obj) <- vln_by
  d <- file.path(base_dir, out_dir, "qc_metrics", "violin_plots")
  dir.create(d, showWarnings = FALSE)
  
  feats <- c("nCount_RNA", "nFeature_RNA")
  for (elem in feats) {
    p <- VlnPlot(seurat_obj, features = elem, pt.size = 0)
    ggsave(p, filename = file.path(d, glue("{elem}.png")), device = "png")
  }
}

# TODO: Finish description

# Example expected structure of ref_genes.csv:
#         ref_gene_colname          ref_cell_colname
#      "PDGFRB, FOXP1, CD34, ..."   "Purkinje Cell"
#       "COL1, COL12, FLT1..."      "Epithelial Cell" ...

# ref_gene_colname: column name containing string of genes corresponding to
#                   cell type of interest
# ref_cell_colname: column name containing cell names as defined by the genes
#                   in ref_gene_colname
# obj:              seurat object of interest, on which FindMarkers has been run

#' findmarkers_gene2cell_mapping
#' @import Seurat glue
#' @export
findmarkers_gene2cell_mapping <- function(obj, ref_genes_csv, findmarkers_csv,
                                          resolutions = seq(0.1, 1, by = 0.1), ref_gene_colname, ref_cell_colname,
                                          out_dir = glue("{getwd()}/FindMarkers_gene2cell_mapping"),
                                          assay_snn = "RNA_snn_res") {
  ref_genes <- read.csv(ref_genes_csv)
  cell_marks <- ref_genes[[ref_gene_colname]]
  cell_names <- ref_genes[[ref_cell_colname]]
  for (r in resolutions){
    Idents(obj) <- glue("{assay_snn}.{r}")
    for (clust in levels(obj)) {
      # load in findmarkers_run
      findmarkers_genes <- read.csv(findmarkers_csv)
      findmarkers_genes <- findmarkers_genes[order(findmarkers_genes$avg_log2FC,
                                                   decreasing = TRUE), ]
      findmarkers_genes <- findmarkers_genes$X[1:30]
      # initialize
      important_gene <- c()
      gene_index <- c()
      cell_type_index <- c()
      i <- 1
      for (elem in cell_marks) {
        genes <- strsplit(x = elem, split = ", ")
        for (g in genes[[1]]) {
          if (g %in% findmarkers_genes) {
            important_gene[length(important_gene) + 1] <- g
            gene_index[length(gene_index) + 1] <- which(g == findmarkers_genes)
            cell_type_index[length(cell_type_index) + 1] <- cell_names[i]
          }
        }
        i <- i + 1
      }
      # if there were any matches between g and findmarkers_genes...
      if (!is.null(important_gene)) {
        # build dataframe to save as csv
        complete_index <- data.frame(important_gene, gene_index,
                                     cell_type_index)
        colnames(complete_index) <- c("important_findmarker_gene",
                                      "avglog2FC_index", "associated_cell_name")
        d <- file.path(out_dir, glue("res{r}"))
        dir.create(d, recursive = TRUE, showWarnings = FALSE)
        write.csv(complete_index, file.path(d, glue("{clust}.csv")))
      }
    }
  }
}


# TODO: add description. Prob should be moved to utilities file
#       (integration package), also include optional threshold
# if family of genes, supply pattern
make_expression_overlays <- function(obj, genes, group_by, outdir,
                                     pattern = "^", threshold = 0,
                                     verbose = FALSE) {
  selected_cells <- c()
  for (gene in genes_of_interest) {
    p1 <- DimPlot(obj, reduction = "umap", group.by = group_by, pt.size = 0.5) +
          ggtitle(gene) + NoAxes() + NoLegend()
    # Find autoset ranges
    xrange <- layer_scales(p1)$x$range$range
    yrange <- layer_scales(p1)$y$range$range
    print(gene)
    if (startsWith(gene, pattern)) {
      # Use grepl to identify columns (genes) that match the pattern
      matching_genes <- rownames(obj)[grepl(gene, rownames(obj))]
      # print(matching_genes)
      for (matched in matching_genes) {
        before <- length(selected_cells)
        selected_cells <- c(selected_cells,
                  colnames(obj)[obj@assays$RNA@data[matched, ] > threshold])
        selected_cells <- unique(selected_cells) # prune cells already on list
        after <- length(selected_cells)
        if (verbose) {
          message(print(glue("Num cells containing gene {matched}:
                             {after - before}")))
        }
      }
      } else {
          # if gene in rownames
          if (gene %in% rownames(obj)) {
            before <- length(selected_cells)
            # Create a logical condition for gene expression
            selected_cells <-
                        colnames(obj)[obj@assays$RNA@data[gene, ] > threshold]
            selected_cells <- unique(selected_cells) # no duplicates
            after <- length(selected_cells)
            if (verbose) {
              message(print(glue("Num cells containing gene {matched}: {after -
                                 before}")))
            }
            if (length(selected_cells) > 0) {
              # Create a subset of the VLMC_brain object containing only the
              # selected cells
              subset_obj <- obj[ ,selected_cells]
              # plot
              p2 <- DimPlot(subset_obj, reduction = "umap", label = FALSE) +
                        NoLegend() + ggtitle(gene) +
                        theme(plot.title = element_text(hjust = 0.5)) +
                        xlim(xrange) + ylim(yrange) + NoAxes()
              im <- integration::magick_overlay(main_plt = p2, trans_plt = p1,
                                  dest = outdir,
                                  glue("by_{my_slot}"),
                                  filename = glue("{gene}.png"))
            } else {
              message(print(glue("Num cells expressing {gene} does not meet
                                 provided threshold: {threshold}.")))
            }
          } else {
            message(print(glue("{gene} is not a gene name present in Seurat
                               object.")))
          }
      }
  }
}