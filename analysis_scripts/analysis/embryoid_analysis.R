#!/usr/bin/env Rscript

# Purpose:
#   TODOs included...

library(Matrix)
npcs <- 30

# set paths
base_dir <- getwd()
data_dir <- file.path(base_dir, "data", "embryoid") # premade, stores raw data

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

# perform preprocessing
embryoid <- FindVariableFeatures(embryoid)
embryoid <- ScaleData(embryoid)
embryoid <- NormalizeData(embryoid)
embryoid <- RunPCA(embryoid, assay = "RNA", npcs = npcs)
embryoid <- FindNeighbors(embryoid, dims = 1:npcs, assay = "RNA",
                          verbose = FALSE)
embryoid <- RunUMAP(embryoid, dim = 1:npcs)
embryoid <- FindClusters(embryoid)




# below is null because not shipped as RDS...
command <- intersect(x = c("FindIntegrationAnchors", "FindTransferAnchors"),
                     y = Command(object = embryoid))[1]

Command(object = embryoid, command = command, value = "normalization.method")
