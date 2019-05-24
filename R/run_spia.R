#' This script creates pathway files for SPIA
#' from excel files

###~~~~~main function to read all sheets in a excel file~~~~~###
read_allsheets <- function(filename, tibble = FALSE) {
  require(readxl)
  sheets <- readxl::excel_sheets(filename)
  tmp <- lapply(sheets, function(x) readxl::read_excel(filename, sheet = x,col_types ="numeric",col_names = TRUE))
  if(!tibble) tmp <- lapply(tmp, as.matrix)
  tmp <- lapply(tmp,function(x) {rownames(x) <- colnames(x); x})
  names(tmp) <- sheets
  tmp
}

##~~~read all excel fles in a folder to create te main lsit objects~~~~~##
##~~~keep all excel files in a single folder~~~~~~##
##~~~ change the directory accordingly each time 

setwd("/Users/danieldomingo/.pathme/spia/kegg")
#setwd(/Users/danieldomingo/.pathme/spia/reactome")
#setwd("/Users/danieldomingo/.pathme/spia/wikipathways")

# Read the files and munges every identified from the name
filenames <-list.files(path = ".")
path.info = lapply(paste0(getwd(),"/", filenames), read_allsheets) 
filenames = gsub('_unflatten','',filenames)
filenames = gsub('.xlsx','',filenames)

names(path.info) = filenames

sprintf("Total number of pathways to be analyzed: %s", length(filenames))


###~~~~~~add last 3 objects in each list inside the path.info~~~~###

# Prepare SPIA data and stores it as a RData file
for (pathway in 1:length(path.info)) {
  for (matrix in 1:length(path.info[pathway])) {
    # Required data structure by SPIA
    path.info[[pathway]][['nodes']] <- rownames(path.info[[pathway]][[matrix]])
    path.info[[pathway]][['title']] <- filenames[pathway]
    path.info[[pathway]][['NumberOfReactions']] <- as.integer(0)
    names(path.info[[pathway]])[[3]] <- "binding/association"
    names(path.info[[pathway]])[[22]] <- "activation_binding/association"
  }
}

##############################################
save(path.info,file="/Users/danieldomingo/Downloads/hsaSPIA.RData")
#####################


#####________SPIA____________#####
# Load required libraries
library(SPIA)
library(tidyr)

#' Read the differential expression data
data = read.csv('/Users/danieldomingo/Downloads/prad_deseq2.csv')
#' Extract two columns we need (HGNC Symbol and log2foldchange)
data = data[,c(8,3)]

#' Read the gene universe file (list of HGNC symbols)
universe = read.csv('/Users/danieldomingo/PycharmProjects/compath-revolutions/spia/hgnc_universe.txt', header = 0,stringsAsFactors = F)
names(universe) = "gene_symbol"

#' merge data with the universe to get only those genes which are in universe
tmp = merge(data, universe,by="gene_symbol")
tmp = unique(tmp)
# Make sure there are no NAs
tmp = tmp %>% drop_na()

#' create a named vector of the dataset and for the universe with the right format (character and numeric)
myData = setNames(as.numeric(tmp$log2FoldChange),as.character(tmp$gene_symbol))
myUniverse = as.character(universe$gene_symbol)

spia_wiki = spia(
  de = myData,
  all = myUniverse,
  data.dir="/Users/danieldomingo/Downloads/", organism="hsa"
)

setwd("/Users/danieldomingo/Downloads")
# Write CSV in R
write.csv(spia_wiki, file = "prad_reactome_spia.csv")

#names(spia_wiki) = names(spia_wiki)
#save(spia_wiki,file = "/home/memon/projects/msdrp/data/spia/spia_wikipath_custom_results.RData")
# save(spia_kegg,file = "/home/memon/projects/msdrp/data/spia/spia_kegg_custom_results.RData")

