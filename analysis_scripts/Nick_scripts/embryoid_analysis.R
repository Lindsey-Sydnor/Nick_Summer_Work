# AUTHOR: Nick Sibiryakov

# important libraries
library(glue)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(clustree)
library(ggraph)
library(Matrix)
library(stringr)
# set a base directory
base_dir <- "/home/nsibir/oxford_organoid_testing"
dirs <- c("images", "objs")
for (d in dirs){
    dir.create(file.path(base_dir, d), recursive = TRUE, showWarnings = FALSE)
}
image_dir <- file.path(base_dir, "images")
obj_dir <- file.path(base_dir, "objs")
umap_dir <- file.path(base_dir, "umaps")
# analyze files
fullmatrix <- readMM(file.path(glue("/home/nsibir/oxford_organoid_testing/",
                                            "objs/GSM3573649_D_matrix.mtx")))
barcodes <- read.table(file = file.path(glue(base_dir,
                                         "/objs/GSM3573649_D_barcodes.tsv")))
genes <- read.table(file = file.path(glue(base_dir,
                                            "/objs/GSM3573649_D_genes.tsv")))
# make the matrix actually readable
barcodes <- barcodes$V1
gene_names <- genes$V2
colnames(fullmatrix) <- barcodes
rownames(fullmatrix) <- gene_names
# load in the embryoid data set

embryonic <- CreateSeuratObject(counts = fullmatrix, project = "embryonic4",
                                        min.cells = 3, min.features = 200)
# update seurat object because it's kinda old
embryonic <- UpdateSeuratObject(embryonic)
#find variable features
embryonic <- FindVariableFeatures(embryonic)
#scale data
embryonic <- ScaleData(embryonic)
# normalize data
embryonic <- NormalizeData(embryonic)
# set ncps to 30
npcs <- 30
# run pca 
embryonic <- RunPCA(embryonic, assay = "RNA", npcs = npcs)
# find the neighbors
embryonic <- FindNeighbors(embryonic, dims = 1:npcs, assay = "RNA", 
                            verbose = FALSE)
# run umap on first five pcs? 
embryonic <- RunUMAP(embryonic, dim = 1:5)
# find clusters
embryonic <- FindClusters(embryonic)
# add to meta data
cell_neighbors <- Idents(embryonic)

# calculate mitochondrial dna
mito.genes = grep(pattern = "^MT-", x = rownames(embryonic), value = TRUE)
percent.mito = Matrix::colSums(embryonic[mito.genes, 
                                                  ])/Matrix::colSums(embryonic)
# add mitochondrial dna to meta data
embryonic = AddMetaData(object = embryonic, metadata = percent.mito, 
                                                    col.name = "percent.mito")
# feature plot
p <- FeaturePlot(embryonic, features = "percent.mito")
ggsave(filename = file.path(base_dir, "images", "mito.png"), plot = p)
# dim plot
p <- DimPlot(embryonic, reduction = "umap")
ggsave(filename = file.path(base_dir, "images", "dim.png"), plot = p)
# overlay mitochondrial onto umap
library(magrittr)
data1 = CreateSeuratObject(counts = embryonic[["RNA"]]@counts[,-7:-19],project="mito")
data2 = CreateSeuratObject(counts = embryonic[["RNA"]]@counts[,41:80],project="ptx")

#violin plot!
p <- VlnPlot(embryonic, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3)
ggsave(filename = file.path(base_dir, "images", "image1.png"), plot = p)

# reading in the cell markers we're looking for
read_csv <- read.csv(glue("/home/nsibir/oxford_organoid_testing/", 
                    "gene_expression_matrix2.csv")) 
cell_marks <- read_csv$X.1[3:length(read_csv$X.1)]
cell_names <- read_csv$X[3:length(read_csv$X)]
#save updated object
saveRDS(object = embryonic, file = file.path(obj_dir, "updated_embryoid.rds"))

# rerun processing pipeline @ multiple resolutions on integrated assay:
resolutions <- seq(0.1, 1, by = 0.1)

for (r in resolutions) {
    embryonic <- FindClusters(embryonic, verbose = FALSE,
                                       resolution = r)
    p <- DimPlot(embryonic, reduction = "umap",
                 group.by = glue("RNA_snn_res.{r}")) + NoAxes()
    ggsave(filename = file.path(umap_dir, glue("embryonic {r}.png")), plot = p)
}
# run find markers on 0.4 resolution and 0.5 to see
# Run FindMarkers on each cluster of each resolution. Save to csv.
for (r in resolutions) {
  # set active ident to r of interest
  Idents(embryonic) <- glue("RNA_snn_res.{r}")
  print(glue("Running resolution {r}..."))

  # make a directory for r (each resolution)
  d <- file.path(obj_dir, "FindMarkers_run", glue("res{r}"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

  # FindMarkers for each cluster
  for (clust in (levels(embryonic))){
    print(glue("Running FindMarkers on cluster {clust}..."))
    cell_markers <- FindMarkers(embryonic, ident.1 = clust,
                                test.use = "bimod", min.pct = 0.1,
                                only.pos = TRUE,
                                max.cells.per.ident = Inf)
    cellmark_df <- as.data.frame(cell_markers)
    # save to csv
    write.csv(x = cellmark_df, file = file.path(d, glue("{clust}.csv")))
  }
}
##### Comparing externally-defined cell type markers with FindMarkers output


for (r in resolutions) {
  Idents(embryonic) <- glue("RNA_snn_res.{r}")
  for (clust in levels(embryonic)) {
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

# run clustree to get a sense of the ideal resolution
library(clustree)
library(ggraph)
library(ggplot2)
clust_diagram <- clustree(embryonic, prefix = "RNA_snn_res.")
ggsave(filename = file.path(image_dir, "clustree_emb.png"), 
                                                        plot = clust_diagram)

# make a dir to hold res-level dot plots
res_dir <- file.path(image_dir, "embryoid_dotplots")
dir.create(res_dir)
# create res-level dot plots
for (r in resolutions) {
  important_gene <- c()
  for (elem in cell_marks) {
    gene <- strsplit(elem, split = ", ")
    for (g in gene[[1]]) {
      important_gene[length(important_gene) + 1] <- g
    }
  }
  # to see if any important genes weren't in our data set
  any(!(important_gene %in% rownames(embryonic)))

  # keep only unique values
  important_gene <- unique(important_gene)

  p <- DotPlot(embryonic, assay = "RNA", features = important_gene,
                  group.by = glue("RNA_snn_res.{r}")) +
          labs(title = str_wrap(glue("Resolution {r}"), 35)) +
          xlab("Genes") +
          ylab("Cluster Identity") +
          theme(axis.text.x = element_text(angle = 90, size = 7))

  ggsave(plot = p, filename = file.path(res_dir, glue("{r}.png")))
}

# Initialize an empty data frame with column names
gene_markers <- read.csv("/home/nsibir/oxford_organoid_testing/cbl.54.csv")
results <- matrix(nrow = length(gene_markers$CBL), ncol = length(unique(embryonic$cell_classes)))
rownames(results) <- gene_markers$CBL
colnames(results) <- unique(aldinger_organoid2$cell_classes)

options(scipen = 100, digits = 4) # set global digit num


results_dir <- "/home/nsibir/oxford_organoid_testing"
for (cell in colnames(results)) {
  # if genes in average_expression are in gene_markers$CBL 
#   subsets the genes in average_expression
  matched <- average_expression[which(rownames(average_expression) %in% gene_markers$CBL), cell]
  # matches the names of the genes in "gene_markers$CBL" to the matched genes
  matched_indices <- match(gene_markers$CBL, names(matched))
  # per cell column it populates all gene expression based off matching indices with gene_marker$CBL
  # to make this in the order of gene_marker$CBL 
  results[, cell] <-  as.numeric(matched[matched_indices])
}
 write.csv(results, file.path(results_dir, "gene_expression_matrix2.csv"),
                row.names = TRUE)














###########################################################################

# RunPCA & create nearest neighbor (NN) / shared NN graphs for RNA assay
ncps <- 30
embryonic <- RunPCA(embryonic, assay = "RNA", npcs = npcs)
embryonic <- FindNeighbors(embryonic, dims = 1:npcs, assay = "RNA", 
                            verbose = FALSE)
# find average expression
average_expression <- AverageExpression(object = embryonic,
                                        assays = "RNA")        
average_expression <- average_expression[[1]]      

gene_markers <- read.csv("/home/nsibir/gene_expression/gene.data.csv")
results <- matrix(nrow = length(gene_markers$CBL)
results_dir <- "/home/nsibir/oxford_organoid_testing"
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

results_dir <- "/home/nsibir/oxford_organoid_testing"
for (cell in colnames(results)) {
  # if genes in average_expression are in gene_markers$CBL subsets the genes in average_expression
  matched <- average_expression[which(rownames(average_expression) %in% gene_markers$CBL), cell]
  # matches the names of the genes in "gene_markers$CBL" to the matched genes
  matched_indices <- match(gene_markers$CBL, names(matched))
  # per cell column it populates all gene expression based off matching indices with gene_marker$CBL
  # to make this in the order of gene_marker$CBL 
  results[, cell] <-  as.numeric(matched[matched_indices])
}
 write.csv(results, file.path(results_dir, "embryonic_gene_expression_matrix.csv"),
                row.names = TRUE)

                                
# average expression for each cell type
average_expression <- AverageExpression(object = embryonic,
                                        group.by = "ident",
                                        assays = "RNA")
average_expression <- average_expression[[1]]
# enrichment analysis
library(escape)
library(Seurat)
library(dittoSeq)
library(ggplot2)
library(GSEABase)
library(GSVA)
# color stuff
colorblind_vector <- colorRampPalette(rev(c("#0D0887FF", "#47039FFF", 
              "#7301A8FF", "#9C179EFF", "#BD3786FF", "#D8576BFF",
              "#ED7953FF", "#FA9E3BFF", "#FDC926FF", "#F0F921FF")))
# dim plot
c <- DimPlot(embryonic, label = T, reduction = "umap") + NoLegend()
ggsave(plot = c, filename = file.path(image_dir, "embdimplot.png"))
# initializing gene set
gene.sets <- list(important_genes = important_gene)
# run enrichment
ES <- enrichIt(obj = embryonic,
               gene.sets = gene.sets,
               groups = 1000, cores = 4)
#restore meta data
embryonic <- AddMetaData(embryonic, ES)
# enrichment in a ridgeplot 
ES2 <- data.frame(embryonic[[]], Idents(embryonic))
colnames(ES2)[ncol(ES2)] <- "cluster"
p <- ridgeEnrichment(ES2, group = "cluster", gene.set = "important_gene", add.rug = TRUE)
ggsave(plot = p, filename = file.path(image_dir, "enrichplot.png"))
# enrichment plot