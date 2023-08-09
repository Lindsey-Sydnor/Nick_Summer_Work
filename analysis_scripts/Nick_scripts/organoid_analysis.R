# AUTHOR: Nick Sibiryakov, some edits/additions by Lindsey Sydnor

library(Seurat)
library(dplyr)
library(patchwork)
library(glue)
library(ggplot2)
library(clustree)
library(ggraph)

# set base dir:
base_dir <- "/home/nsibir/oxford_organoid_testing/repo.stuff"

# creating output directories
dirs <- c("images", "objs")
for (d in dirs){
  dir.create(file.path(base_dir, d), recursive = TRUE, showWarnings = FALSE)
}

# setting defaults
image_dir <- file.path(base_dir, "images")
obj_dir <- file.path(base_dir, "objs")
umap_dir <- file.path(image_dir, "umaps")

dir.create(umap_dir)

# presets
npcs <- 30

### Run first time:
# read in file with readRDS
aldinger_organoid2 <- readRDS(glue("/home/nsibir/oxford_organoid_testing/",
                              "repo.stuff/data/neworganoids_combined.rds"))

# optional: set active assay to RNA
DefaultAssay(aldinger_organoid2) <- "RNA"

# plot original UMAP
p <- DimPlot(aldinger_organoid2, reduction = "umap")
ggsave(filename = file.path(image_dir, "orig2_umap.png"), plot = p)

#find variable features
aldinger_organoid2 <- FindVariableFeatures(aldinger_organoid2)

#scale data
aldinger_organoid2 <- ScaleData(aldinger_organoid2)

# normalize data
aldinger_organoid2 <- NormalizeData(aldinger_organoid2)

# save weird extra identity, creates a metadata column
cell_classes <- Idents(aldinger_organoid2)
aldinger_organoid2$cell_classes <- gsub(" ", "_", cell_classes)
saveRDS(aldinger_organoid2, glue("/home/nsibir/oxford_organoid_testing/",
                                   "organoids_cell_classes.rds"))
aldinger_organoid2 <- readRDS(glue("/home/nsibir/oxford_organoid_testing/",
                                   "organoids_cell_classes.rds"))

##### Remove mitochondrial and ribosomal genes ##### How to find?
patterns <- c("^MT-", "^MRPL", "^MRPS", "^RPS", "^RPL")
remov_genes <- character()
for (p in patterns) {
  remov_genes <- append(remov_genes, grep(pattern = p,
                        x = rownames(aldinger_organoid2),
                        value = TRUE), length(remov_genes))
}

# save to your altered data folder
saveRDS(object = aldinger_organoid2, file = file.path(obj_dir,
                                                      "updated_organoid.rds"))

############# checkpoint 1 #############
aldinger_organoid2 <- readRDS(file.path(obj_dir, "updated_organoid.rds"))

#create a dimplot and group by "cell classes"
p <- DimPlot(aldinger_organoid2, reduction = "umap", group.by = "cell_classes")
ggsave(filename = file.path(image_dir, "cell_classes_umap1.png"), plot = p)

# check their default cluster labels (not cell type)
length(unique(aldinger_organoid2$seurat_clusters)) ==
    length(unique(aldinger_organoid2$cell_classes)) # TRUE

# check to see if counts are normalized (looking for floats)
aldinger_organoid2@assays$RNA@data # floats

#elbowplot
b <- ElbowPlot(aldinger_organoid2)
ggsave(filename = file.path(base_dir, "images", "elbowplot.png"), plot = b)

#violin plot of mitotic dna ?? for fun?
c <- VlnPlot(aldinger_organoid2, features = 'percent.mt')
ggsave(filename = file.path(base_dir, "images", "percentmt.png"), plot = c)

# rerun processing pipeline @ multiple resolutions on integrated assay:
resolutions <- seq(0.1, 1, by = 0.1)

# RunPCA & create nearest neighbor (NN) / shared NN graphs for RNA assay
aldinger_organoid2 <- RunPCA(aldinger_organoid2, assay = "RNA", npcs = npcs)
aldinger_organoid2 <- FindNeighbors(aldinger_organoid2, dims = 1:npcs,
                                    assay = "RNA", verbose = FALSE)

# Find clusters @ each resolution
for (r in resolutions) {
    aldinger_organoid2 <- FindClusters(aldinger_organoid2, verbose = FALSE,
                                       resolution = r)
    p <- DimPlot(aldinger_organoid2, reduction = "umap",
                 group.by = glue("RNA_snn_res.{r}")) + NoAxes()
    ggsave(filename = file.path(umap_dir, glue("{r}.png")), plot = p)
}
head(pbmc_small[["RNA"]][[]])

############# checkpoint 2 #############
saveRDS(object = aldinger_organoid2,
        file = file.path(obj_dir, "updated_organoid2.rds"))
aldinger_organoid2 <- readRDS(file.path(obj_dir, "updated_organoid2.rds"))
# run clustree for ideal resolution
clust_diagram <- clustree(aldinger_organoid2, prefix = "RNA_snn_res.")
ggsave(filename = file.path(image_dir, "clustree.png"), plot = clust_diagram)

# Run FindMarkers on each cluster of each resolution. Save to csv.
for (r in resolutions) {
  # set active ident to r of interest
  Idents(aldinger_organoid2) <- glue("RNA_snn_res.{r}")
  print(glue("Running resolution {r}..."))

  # make a directory for r (each resolution)
  d <- file.path(obj_dir, "FindMarkers_run", glue("res{r}"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

  # FindMarkers for each cluster
  for (clust in (levels(aldinger_organoid2))){
    print(glue("Running FindMarkers on cluster {clust}..."))
    cell_markers <- FindMarkers(aldinger_organoid2, ident.1 = clust,
                                test.use = "bimod", min.pct = 0.1,
                                only.pos = FALSE,
                                max.cells.per.ident = Inf)
    cellmark_df <- as.data.frame(cell_markers)
    # save to csv
    write.csv(x = cellmark_df, file = file.path(d, glue("{clust}.csv")))
  }
}

##### Comparing externally-defined cell type markers with FindMarkers output

# read in csv file
read_csv <- read.csv(glue("/home/nsibir/oxford_organoid_testing/objs/",
                          "cellmarkers.csv"))
cell_marks <- read_csv$X.1[3:length(read_csv$X.1)]
cell_names <- read_csv$X[3:length(read_csv$X)]

for (r in resolutions) {
  Idents(aldinger_organoid2) <- glue("RNA_snn_res.{r}")
  for (clust in levels(aldinger_organoid2)) {
    # load in findmarkers_run
    findmarkers <- read.csv(glue("/home/nsibir/oxford_organoid_testing/objs/",
                                "FindMarkers_run/res{r}/{clust}.csv"))
    findmarkers <- findmarkers[order(findmarkers$avg_log2FC,
                                     decreasing = TRUE), ]
    findmarkers_genes <- findmarkers$X[1:30]
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
      complete_index <- data.frame(important_gene, gene_index, cell_type_index)
      colnames(complete_index) <- c("important_findmarker_gene",
                                  "avglog2FC_index", "associated_cell_name")
      d <- file.path(obj_dir, "FindMarkers_gene2cell_mapping", glue("res{r}"))
      dir.create(d, recursive = TRUE, showWarnings = FALSE)
      write.csv(complete_index, file.path(d, glue("{clust}.csv")),
                row.names = FALSE)
    }
  }
}
# make a dir to hold res-level dot plots
res_dir <- file.path(image_dir, "res_marker_dotplots")
dir.create(res_dir)
for (r in resolutions) {
  important_gene <- c()
  for (elem in cell_marks) {
    gene <- strsplit(elem, split = ", ")
    for (g in gene[[1]]) {
      important_gene[length(important_gene) + 1] <- g
    }
  }
  # to see if any important genes weren't in our data set
  any(!(important_gene %in% rownames(aldinger_organoid2)))

  # keep only unique values
  important_gene <- unique(important_gene)

  p <- DotPlot(aldinger_organoid2, assay = "RNA", features = important_gene,
                  group.by = glue("RNA_snn_res.{r}")) +
          labs(title = str_wrap(glue("Resolution {r}"), 35)) +
          xlab("Genes") +
          ylab("Cluster Identity") +
          theme(axis.text.x = element_text(angle = 90, size = 7))

  ggsave(plot = p, filename = file.path(res_dir, glue("{r}.png")))
}

#average expression for each cell type
results_dir <- "/home/nsibir/oxford_organoid_testing/repo.stuff/objs"
# Initialize an empty data frame with column names
gene_markers <- read.csv("/home/nsibir/oxford_organoid_testing/cbl.54.csv")
aldinger_organoid2$cell_classes <- gsub(" ", "_", aldinger_organoid2$cell_classes)
results <- matrix(nrow = length(gene_markers$CBL), ncol = length(unique(aldinger_organoid2$cell_classes)))
rownames(results) <- gene_markers$CBL
colnames(results) <- unique(aldinger_organoid2$cell_classes)
options(scipen = 100, digits = 4) # set global digit num
average_expression <- AverageExpression(object = aldinger_organoid2,
                                        group.by = "cell_classes",
                                        assays = "RNA")
average_expression <- average_expression[[1]]
# filter out null gene expression
# Initialize an empty list to store results
average_expression_per_class <- list()
# Iterate over cell classes
for (cell_class in colnames(average_expression)) {
  # Filter expression for the current cell class
  expression_subset <- average_expression[, cell_class]
  # Calculate average expression only if gene expression is > 0
  if (any(expression_subset > 0, na.rm = TRUE)) { #mean(expression..., na.rm = TRUE) 
    average_expression_per_class[[cell_class]] <- expression_subset[expression_subset > 0]
  }
}
average_expression_per_class <- average_expression_per_class[[1]]
average_expression_per_class <- as.matrix(average_expression_per_class)
for (cell in colnames(results)) {
  # if genes in average_expression are in gene_markers$CBL subsets the genes in average_expression
  matched <- average_expression_per_class[which(rownames(average_expression_per_class) %in% gene_markers$CBL), cell]
  # matches the names of the genes in "gene_markers$CBL" to the matched genes
  matched_indices <- match(gene_markers$CBL, names(matched))
  # per cell column it populates all gene expression based off matching indices with gene_marker$CBL
  # to make this in the order of gene_marker$CBL 
  results[, cell] <-  as.numeric(matched[matched_indices])
}
 write.csv(results, file.path(results_dir, "gene_expression_matrix2.csv"),
                row.names = TRUE)
colnames(average_expression) <- gsub(" ", "_", colnames(average_expression))

# ADD CBL COLUMN
for (cell in colnames(results)) {
  # if genes in average_expression are in gene_markers$CBL subsets the genes in average_expression
  matched <- average_expression[which(rownames(average_expression) %in% gene_markers$CBL), cell]
  # matches the names of the genes in "gene_markers$CBL" to the matched genes
  matched_indices <- match(gene_markers$CBL, names(matched))
  # per cell column it populates all gene expression based off matching indices with gene_marker$CBL
  # to make this in the order of gene_marker$CBL 
  results[, cell] <-  as.numeric(matched[matched_indices])
}
 write.csv(results, file.path(results_dir, "gene_expression_matrix2.csv"),
                row.names = TRUE)
# keep all genes that are less than 1
for (i in 1:nrow(results)){
  for (j in 1:ncol(results)) {
    if (results[i, j] > 1){
      results[i, j] <- NA
    }
  }
}

# create a heatmap of the organoids in the results matrix!
# omit na values first!
results[is.na(results)] <- 0
library(pheatmap)
g <- pheatmap(results, cluster_cols = TRUE)
ggsave(plot = g, filename = file.path(image_dir, "final_heatmap.png"))

# Create a Seurat object from the 'results' matrix
seurat_obj <- CreateSeuratObject(counts = results)

# Normalize and scale the data
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

# Create the heatmap using DoHeatmap
c <- DoHeatmap(aldinger_organoid2, features = important_gene, 
                  group.by = "cell_classes") + 
                  theme(text = element_text(size = 10))
ggsave(plot = c, filename = file.path(image_dir, "gene_markers_heatmap.png"))

#load data
data <- read.delim(aldinger_organoid2, header=T, row.names="gene")

f <- aldinger_organoid2@assays$RNA@data[rownames(aldinger_organoid2) %in%
                                        gene_markers$CBL, ]
rownames(aldinger_organoid2) %in% gene_markers$CBL

# heatmap using pheatmap
g <- pheatmap(f, )
ggsave(plot = g, filename = file.path(results_dir, "actualheatmapcbl.png"))

subsetted_cbl <- gene_markers$CBL[!sapply(gene_markers$CBL, identical, "TUBA1A")]

subsetted_cbl <- gene_markers$CBL[gene_markers$CBL != "TUBA1A" &
                                  gene_markers$CBL != "TUBB2A" &
                                  gene_markers$CBL != "CDC42" &
                                  gene_markers$CBL != "AP1S2"]

mat <- average_expression[rownames(average_expression) %in% subsetted_cbl, ]

# heatmap using averageexpression or pseudobulk
g <- pheatmap(mat, cluster_cols = TRUE)
ggsave(plot = g, filename = file.path(results_dir, "avgexheatmap2.png"))
#######CHECKPOINT 3###############
saveRDS(object = aldinger_organoid2,
        file = file.path(obj_dir, "updated_organoid3.rds"))
aldinger_organoid2 <- readRDS(file.path(obj_dir, "updated_organoid3.rds"))

##ssGSEA

library(escape)
library(Seurat)
library(dittoSeq)
library(ggplot2)
library(GSEABase)
library(GSVA)

colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
              "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
              "#ED7953FF", "#FA9E3BFF", "#FDC926FF", "#F0F921FF")))

c <- DimPlot(aldinger_organoid2, label = T, reduction = "umap") + NoLegend()
ggsave(plot = c, filename = file.path(results_dir, "2dimplot.png"))

Idents(aldinger_organoid2) <- "cell_classes"

gene.sets <- list(important_genes = important_gene)

# run the enrichment analysis!!
ES <- enrichIt(obj = aldinger_organoid2,
               gene.sets = gene.sets,
               groups = 1000, cores = 4)

#restore meta data
aldinger_organoid2 <- AddMetaData(aldinger_organoid2, ES)

#omit the na
# aldinger_organoid3 <- complete.cases(aldinger_organoid2)
#create a ditto heatmap

cell_names <- colnames(aldinger_organoid2)
cells_classes <- aldinger_organoid2$cell_classes

# ridge plot
ES2 <- data.frame(aldinger_organoid2[[]], Idents(aldinger_organoid2))
colnames(ES2)[ncol(ES2)] <- "cluster"
p <- ridgeEnrichment(ES2, gene.set = "important_genes", group = "cluster",
                     add.rug = TRUE)
ggsave(plot = p, filename = file.path(image_dir, "ridgeplots.png"))

# why not a violin plot too
l <- dittoPlot(aldinger_organoid2, var = "important_genes",
                group.by = "cell_classes") +
    # theme(legend.position = "none") + # if we don't want legend
    scale_fill_manual(values = colorblind_vector(12))
ggsave(plot = l, filename = file.path(image_dir, "genex_vioplot.png"))