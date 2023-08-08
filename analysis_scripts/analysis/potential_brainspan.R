# https://github.dev/deevdevil88/10x_Human_Cerebellum_Organoid/blob/main/Organoid_Brainspan_Pseudobulk_data_comparison.Rmd

brainspan_original <- read.csv(file = "./genes_matrix_csv/expression_matrix.csv", sep=",",header=F)
brainspan_colnames <- read.csv(file="./genes_matrix_csv/columns_metadata.csv", sep=",")
brainspan_rownames <- read.csv(file="./genes_matrix_csv/rows_metadata.csv", sep=",")
brainspan_original <- brainspan_original[,-1]

rownames(brainspan_original) <- brainspan_rownames$ensembl_gene_id
library(plyr)
brainspan_colnames$dev_period <- brainspan_colnames$age
#brainspan_colnames$dev_period <- plyr::revalue(brainspan_colnames$dev_period, c("1 yrs"="Childhood","10 mos"="Infancy","11 yrs" ="Childhood","12 pcw"= "Early_prenatal","13 pcw"="Early_mid_prenatal","13 yrs" = "Adolescence","15 yrs"="Adolescence","16 pcw"="Early_mid_prenatal", "17 pcw"="Early_mid_prenatal","18 yrs"="Adolescence","19 pcw"="Late_mid_prenatal", "19 yrs"="Adolescence", "2 yrs"= "Childhood", "21 pcw"="Late_mid_prenatal", "21 yrs"="Adults","23 yrs" ="Adults", "24 pcw" ="Late_mid_prenatal", "25 pcw"="Late_prenatal", "26 pcw"="Late_prenatal", "3 yrs"="Childhood", "30 yrs"="Adults", "35 pcw"="Late_prenatal","36 yrs" ="Adults","37 pcw"="Late_prenatal", "37 yrs"="Adults", "4 mos"= "Infancy", "4 yrs" ="Childhood", "40 yrs"="Adults","8 pcw"="Early_prenatal", "8 yrs" ="Childhood","9 pcw"= "Early_prenatal" ))

brainspan_colnames$dev_period <- plyr::revalue(brainspan_colnames$dev_period, c("1 yrs"="Infancy","10 mos"="Infancy","11 yrs" ="Childhood","12 pcw"= "Early_prenatal","13 pcw"="Early_prenatal","13 yrs" = "Adolescence","15 yrs"="Adolescence","16 pcw"="Early_prenatal", "17 pcw"="Early_prenatal","18 yrs"="Adolescence","19 pcw"="Late_prenatal", "19 yrs"="Adults", "2 yrs"= "Childhood", "21 pcw"="Late_prenatal", "21 yrs"="Adults","23 yrs" ="Adults", "24 pcw" ="Late_prenatal", "25 pcw"="Late_prenatal", "26 pcw"="Late_prenatal", "3 yrs"="Childhood", "30 yrs"="Adults", "35 pcw"="Late_prenatal","36 yrs" ="Adults","37 pcw"="Late_prenatal", "37 yrs"="Adults", "4 mos"= "Infancy", "4 yrs" ="Childhood", "40 yrs"="Adults","8 pcw"="Early_prenatal", "8 yrs" ="Childhood","9 pcw"= "Early_prenatal" ))




# renames level 1 cell suptyes to major cell classes
brainspan_colnames$brain_region_short <- brainspan_colnames$structure_acronym
brainspan_colnames$brain_region_short <-  plyr::revalue(brainspan_colnames$brain_region_short, c("A1C"="Cortical", "AMY"="Amygdaloid","CB"="Cerebellum","CBC"="Cortical","CGE"="Ganglionic_eminence","DFC"="Cortical","DTH"="Thalamus","HIP"="Hippocampus","IPC"="Cortical","ITC"="Cortical","LGE"="Ganglionic_eminence","M1C"="Cortical","M1C-S1C"="Cortical","MD"="Thalamus","MFC"="Cortical", "MGE"="Ganglionic_eminence","Ocx"="Cortical","OFC"="Cortical","PCx"="Cortical","S1C"="Cortical","STC"="Cortical","STR"="Striatum","TCx"="Cortical", "URL"="URL","V1C"="Cortical","VFC"="Cortical"))

brainspan_colnames$Sample_ID  <- paste(brainspan_colnames$structure_acronym,"_",brainspan_colnames$dev_period,sep="")
write.table(brainspan_colnames, file="brainspan_colnames_extended.txt", sep="\t",quote=F)

colnames(brainspan_original) <- brainspan_colnames$Sample_ID

# First average by Dev period and then by Tissue
#brainspan_original <- read.table(file="./brainspan_original_extended.txt",header=T)
brainspan_original_t <- as.data.frame(t(brainspan_original))
brainspan_original_t$Sample_ID <- brainspan_colnames$Sample_ID
#brainspan_original_t$Sample_ID <- gsub('\\..*', '', brainspan_original_t$Sample_ID, perl=T)
#brainspan_original_t$Sample_ID <- rownames(brainspan_original_t)
brainspan_original_t$Sample_ID <- as.factor(brainspan_original_t$Sample_ID)
brainspan_original_t_1 = brainspan_original_t %>%
                  group_by(Sample_ID) %>%
                  summarise_all(mean) 
brainspan_original_t_1 <- as.data.frame(brainspan_original_t_1)
brainspan_original_t_1 <- brainspan_original_t_1[match(brainspan_metadata$Sample_ID, brainspan_original_t_1$Sample_ID),]
rownames(brainspan_original_t_1) <- brainspan_original_t_1$Sample_ID
#drops <- "Sample_ID"
brainspan_original_t_1 <- subset(brainspan_original_t_1, select=-c(Sample_ID))
rm(brainspan_original_t)

brainspan_original <- as.data.frame(t(brainspan_original_t_1))
write.table(brainspan_original, file="brainspan_processed_extended_prenatal.txt", sep="\t",quote=F)
rm(brainspan_original_t_1)

brainspan_metadata <- brainspan_colnames %>% distinct(Sample_ID, .keep_all = TRUE)
write.table(brainspan_metadata, file="brainspan_metadata_extended_prenatal.txt", sep="\t", quote=F)



# read in processed for comparison
brainspan_processed <- read.table(file="./brainspan_processed", sep="\t",header=T, row.names=1)
brainspan_metadata <- as.data.frame(colnames(brainspan_processed))
library(stringr)
brainspan_metadata_1 <-  str_split_fixed(brainspan_metadata$`colnames(brainspan_processed)`, "_", 2)
brainspan_metadata_1 <- as.data.frame(brainspan_metadata_1)
colnames(brainspan_metadata_1)[1] <- "brain_region"
colnames(brainspan_metadata_1)[2] <- "maturity_stage"

library(plyr)
brainspan_metadata_1$brain_region_short <- brainspan_metadata_1$brain_region
# renames level 1 cell suptyes to major cell classes
brainspan_metadata_1$brain_region_short <-  revalue(brainspan_metadata_1$brain_region_short, c("A1C"="Cortical", "AMY"="Amygdaloid","CB"="Cerebellum","CBC"="Cortical","CGE"="Ganglionic_eminence","DFC"="Cortical","DTH"="Thalamus","HIP"="Hippocampus","IPC"="Cortical","ITC"="Cortical","LGE"="Ganglionic_eminence","M1C"="Cortical","M1C.S1C"="Cortical","MD"="Thalamus","MFC"="Cortical", "MGE"="Ganglionic_eminence","Ocx"="Cortical","OFC"="Cortical","PCx"="Cortical","S1C"="Cortical","STC"="Cortical","STR"="Striatum","TCx"="Cortical", "URL"="URL","V1C"="Cortical","VFC"="Cortical"))

brainspan_metadata_1$Sample_ID <- brainspan_metadata$`colnames(brainspan_processed)`
brainspan_metadata_1$Study <- "Brainspan"
#col_order <- c("Sample_ID","Replicate","Group","Study","brain_region")

#brainspan_metadata_1 <- brainspan_metadata_1[,col_order]

rownames(brainspan_metadata_1) <- brainspan_metadata_1$Sample_ID

# read in another person's object
Sams_obj <- readRDS(file="./20190807vs2_v3update_newganoidscombined.rds")
Sams_obj[["cell_type"]] <- Idents(object = Sams_obj)
Sams_obj$cell_type_stim <- paste(Sams_obj$cell_type,"_",Sams_obj$stim,sep="")
Idents(object = Sams_obj) <- "cell_type_stim"
Sams_avg_exp <- AverageExpression(Sams_obj, assay="SCT")
Sams_avg_exp_sct <- as.data.frame(Sams_avg_exp$SCT)
#Sams_avg_exp_sct <- log2(Sams_avg_exp_sct+1)

Sams_avg_exp_sct$gene <- rownames(Sams_avg_exp_sct)
# convert gene symbol to Ens for single cell average expression data
list_gene <- as.character(Sams_avg_exp_sct$gene)
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
query_ensembl<-getBM(attributes=c('ensembl_gene_id','hgnc_symbol','gene_biotype'), mart = ensembl, values=list_gene)

query_ensembl_pc <- subset(query_ensembl, query_ensembl$gene_biotype=="protein_coding", select = 1:2)

Sams_avg_exp_sct_ens <- merge(Sams_avg_exp_sct,query_ensembl_pc,by.x="gene",by.y="hgnc_symbol")
rownames(Sams_avg_exp_sct_ens) <- Sams_avg_exp_sct_ens$ensembl_gene_id
Sams_avg_exp_sct_ens <- subset(Sams_avg_exp_sct_ens, select=-c(gene,ensembl_gene_id))
rm(Sams_avg_exp_sct)





# Remove low expressed genes and only keep protein coding genes in the brainspan
# RPKM matrix
#brainspan_original <- read.table(file="./brainspan_processed_extended_prenatal.txt")
# Remove low expressed genes from brain span RPKM data
list_gene <- as.character(rownames(brainspan_original))
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
query_ensembl<-getBM(attributes=c('ensembl_gene_id','hgnc_symbol','gene_biotype'), mart = ensembl, values=list_gene)

query_ensembl_pc <- subset(query_ensembl, query_ensembl$gene_biotype=="protein_coding", select = 1:3)
brainspan_original$gene <- rownames(brainspan_original)
brainspan_original_pc <- merge(brainspan_original,query_ensembl_pc,by.x="gene",by.y="ensembl_gene_id")
brainspan_original_pc <- subset(brainspan_original_pc, select=-c(hgnc_symbol,gene_biotype))

brainspan_original_pc = brainspan_original_pc %>%
                 group_by(gene) %>%
                summarise_all(mean) 
brainspan_original_pc <- as.data.frame(brainspan_original_pc)

rownames(brainspan_original_pc) <- brainspan_original_pc$gene


brainspan_original_pc <- subset(brainspan_original_pc, select=-c(gene))


#cat("Number of genes and samples: ", dim(brainspan_original_pc), fill=T)
#cat("Number of genes with no read (estimated TPM) in any sample: ", sum(rowSums(brainspan_original_pc) == 0), fill=T)
#cat("Number of genes with exactly a count of 1 in a single sample (genes with row sum of 1):", sum(rowSums(brainspan_original_pc) == 1), fill=T)
# Filtering all genes with o counts across all samples
#keep_feature <- rowSums(brainspan_original_pc > 0) > 0
#brainspan_original_pc <- brainspan_original_pc[keep_feature, ]

# at least 4 samples with a TPM 1 or higher 
#keep_feature <- rowSums(brainspan_original_pc >= 1) >= 3
#brainspan_original_pc <- brainspan_original_pc[keep_feature, ]

#cat("Number of genes and samples after filtering of genes based on a RPKM >= 1 in atleast 3 samples : ", dim(brainspan_original_pc), fill=T)

#brainspan_original_pc <- log2(brainspan_original_pc+1)




# Keep only the common genes between the brainspan and the single cell organoid
# data

#brainspan_metadata_short <- brainspan_metadata_1[,c(4,3,2)]
brainspan_metadata_short <- brainspan_metadata[,c(11,10,9)]
setnames(brainspan_metadata_short,c("brain_region_short","dev_period"), c("cell_type","condition"))

#setnames(brainspan_metadata_short,c("brain_region_short","maturity_stage"), c("cell_type","condition"))
brainspan_metadata_short$Study <- "brainspan"
rownames(brainspan_metadata_short) <- brainspan_metadata_short$Sample_ID
#brainspan_processed <- log2(brainspan_processed+1)



common_genes <- intersect(rownames(brainspan_original_pc), rownames(Sams_avg_exp_sct_ens))

# filter brain span data to common genes
#brainspan_processed$gene <- rownames(brainspan_processed)
#brainspan_processed_common <- brainspan_processed %>% filter(gene %in% common_genes)
#rownames(brainspan_processed_common) <- brainspan_processed_common$gene
#brainspan_processed_common <- subset(brainspan_processed_common, select=-c(gene))

brainspan_original_pc$gene <- rownames(brainspan_original_pc)
brainspan_original_common <- brainspan_original_pc %>% filter(gene %in% common_genes)
rownames(brainspan_original_common) <- brainspan_original_common$gene
brainspan_original_common <- subset(brainspan_original_common, select=-c(gene))



# Filter single cell average expresison object to only the common genes
Sams_avg_exp_sct_ens$gene <- rownames(Sams_avg_exp_sct_ens)
Sams_avg_exp_sct_ens_common <- Sams_avg_exp_sct_ens %>% filter(gene %in% common_genes)
rm(Sams_avg_exp_sct_ens)
rownames(Sams_avg_exp_sct_ens_common) <- Sams_avg_exp_sct_ens_common$gene
Sams_avg_exp_sct_ens_common <- subset(Sams_avg_exp_sct_ens_common, select=-c(gene))

# Create metadata for the single cell Averaged object
Sams_avg_exp <- AverageExpression(Sams_obj, assay="SCT", return.seurat = T)
Sams_avg_exp_metadata <- Sams_avg_exp@meta.data
rm(Sams_avg_exp)
Sams_avg_exp_metadata$stim <- rownames(Sams_avg_exp_metadata)
Sams_avg_exp_metadata <- Sams_avg_exp_metadata %>% separate(stim, sep="_",c("A","B"))
setnames(Sams_avg_exp_metadata, c("A","B"), c("cell_type","condition"))
Sams_avg_exp_metadata$Sample_ID <- rownames(Sams_avg_exp_metadata)

Sams_avg_exp_metadata_short <- Sams_avg_exp_metadata[,c(6,4,5)]
Sams_avg_exp_metadata_short$Study <- "Nayler_Organoid"




# combine Brainspan and Nayler Average expression
# Remove batch effect ysing ComBat

combined_data <- cbind(brainspan_original_common,Sams_avg_exp_sct_ens_common)
combined_data <- log2(combined_data +1)
combined_metadata <- rbind(brainspan_metadata_short,Sams_avg_exp_metadata_short)
rownames(combined_metadata) <- combined_metadata$Sample_ID
combined_metadata$Study <- as.factor(combined_metadata$Study)

library(RColorBrewer)
colB1 <- Vector2Colour(combined_metadata$Sample_ID)
colB2 <- Vector2Colour(combined_metadata$cell_type)
colB3 <- Vector2Colour(combined_metadata$Study)
#colB4 <- Vector2Colour(combined_metadata$Study)

coloresB <- c("colB1","colB2","colB3")
nomsB <- c( "Sample_ID","cell_type","Study")


combat_combined_data <- ComBat(as.matrix(combined_data), batch= combined_metadata$Study, par.prior = T, prior.plots = F)
combat_brainspan <- as.data.frame(combat_combined_data[,1:107])
combat_Sam <- as.data.frame(combat_combined_data[,108:131])

# Run Limma
limma_combined_data <- limma::removeBatchEffect(as.matrix(combined_data), batch=combined_metadata$Study)
limma_brainspan <- limma_combined_data[,1:107]
limma_Sam <- limma_combined_data[,108:131]




# etc.