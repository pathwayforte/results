#' Query and download TCGA datasets from GDC
# R version 3.5.2
# TCGAbiolinks 2.10.3

rm(list = ls());

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

# TODO: Add path to the raw data of the cancer dataset
# Query TCGA-BRCA gene expression quantification from harmonized database 
query <- GDCquery(project = "TCGA-BRCA",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification", 
                       workflow.type = "HTSeq - FPKM-UQ")

# Download and prepare TCGA-BRCA data
GDCdownload(query)
data_expr <- GDCprepare(query)

save(file = "BRCA_RNA-Seq_Data.rda", data_expr)
