#' Preprocess TCGA Data in a way that can be input into PathwayForte

rm(list = ls());

#setwd("/Users/bioadmin/Projects/compath-revolutions/data/tcga_datasets_brca")
install.packages("xlsx", dependencies = TRUE)
library(TCGAWorkflowData) 
library(SummarizedExperiment)
library(RTCGAToolbox)
library(BiocGenerics)
library(Biobase)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(TCGAbiolinks)
library(plyr)
library(dplyr)
library(DT)
library(stringr)
library(xlsx)

############## Prepare expression data as input into GSEA ################

# Input TCGA expression data
# TODO: Add path to the raw data of the cancer dataset
load(file = "./data/BRCA_RNA-Seq_Data.rda")

# TODO: Change the name depending on the loaded dataset
# assay(data) to store expression data
expr_matrix <- assay(brca_data_expr)
#expr_matrix_subset = expr_matrix[1:4,1:2]

# TODO: Change the name depending on the loaded dataset
# Show object data and its metadata
meta_data <- as.data.frame(colData(brca_data_expr))
#transform gene names in data set to HGNC symbols
genes <- rowRanges(brca_data_expr)
gene_id_mappings <- genes$external_gene_name
gene_id_mappings = str_remove(string = gene_id_mappings,pattern = "  ")
rownames(expr_matrix) = genes$external_gene_name

# Get expression values for each P 
expr_data = expr_matrix[gene_id_mappings,]

# Get tissue type from definition column (tumor or normal)
tissue_type = meta_data$definition
tumor = tissue_type == "Primary solid Tumor"
normal = tissue_type == "Solid Tissue Normal"

# Get vital status from vital_status column (alive or dead)
status = meta_data$vital_status

# Get vital status from vital_status column (alive or dead)
days_death = meta_data$days_to_death

# Get matrix tumor, survival and days to death information
merged_survival_info = cbind(tumor, status, days_death)

# Get matrix with columns for survival and status and days to death for tumor cases
event_time_matrix = merged_survival_info[merged_survival_info[, "tumor"] == TRUE,]
event_time_matrix = as.data.frame(event_time_matrix)
event_time_matrix$tumor <- NULL

# Get matrix for normal and tumor gene expression values 
normal_expr = expr_data[,normal]
tumor_expr = expr_data[,tumor]

# Export matrices 
write.table(event_time_matrix, file="event_time_matrix.txt",sep = "\t", row.names=TRUE, col.names=TRUE)
write.table(tumor_expr, file="tumor_expression_matrix.txt",sep = "\t", row.names=TRUE, col.names=TRUE)

#tumor and days_to_death

dim_normal = (dim(normal_expr)[2])
dim_tumor = (dim(tumor_expr)[2])
write.table(dim_normal, file='normal_expression_dimension.txt', sep = "\t", row.names=FALSE, col.names=FALSE)
write.table(dim_tumor, file='tumor_expression_dimension.txt', sep = "\t", row.names=FALSE, col.names=FALSE)

# Join normal and tumor gene expression data matrices
joined_expr_matrix = cbind(normal_expr,tumor_expr)

# TODO: Change the name depending on the dataset in order to export it
write.table(joined_expr_matrix, file="brca_expression_matrix.txt",sep = "\t", row.names=TRUE, col.names=TRUE)
