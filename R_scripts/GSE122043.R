### DA-seq on mouse skin data
### Original paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6361530/
### This script reproduces analysis presented in Figure 3
### https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122043

library(tidyverse)
library(Matrix)
library(reshape2)
library(tidyseurat)
library(Seurat)

setwd("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE122043_mouse_skin/")

##=============================================##
## Functions
# Read in 10X data in .mtx format and add col and row names
read10X <- function(folder = NULL, matFile = NULL, geneFile = NULL, cellFile = NULL, 
                    suffix = "", sep = "_", gz = T){
    if(!is.null(folder)){
        if(gz){
            matFile <- paste(folder, "/matrix.mtx.gz", sep = "")
            geneFile <- paste(folder, "/genes.tsv.gz", sep = "")
            cellFile <- paste(folder, "/barcodes.tsv.gz", sep = "")
        } else {
            matFile <- paste(folder, "/matrix.mtx", sep = "")
            geneFile <- paste(folder, "/genes.tsv", sep = "")
            cellFile <- paste(folder, "/barcodes.tsv", sep = "")
        }
    }
    
    geneNames <- read.table(geneFile, header = F, sep = "\t", as.is = T)[,2]
    cellNames <- paste(read.table(cellFile, header = F, sep = "\t", as.is = T)[,1], suffix, sep = sep)
    
    # add suffix to duplicate gene names
    if(max(table(geneNames)) > 1){
        for(dupgene in names(which(table(geneNames) != 1))){
            geneidx <- which(geneNames == dupgene)
            for(ii in 2:length(geneidx)){
                geneNames[geneidx[ii]] <- paste(dupgene, ii-1, sep = ".")
            }
        }
    }
    
    rawMat <- readMM(matFile)
    rownames(rawMat) <- geneNames
    colnames(rawMat) <- cellNames
    
    return(rawMat)
}

##=============================================##
## Data prep

## Load data
system(
  "tar -xvf ./data/GSE122043_RAW.tar -C ./data/"
)

# Read data of E13.5 and E14.5 samples
sample_names <- c(
  "GSM3453535_e13Control","GSM3453536_e13Control_replicate",
  "GSM3453537_e14Control","GSM3453538_e14Control_replicate",
  "GSM3453539_e13LOF","GSM3453540_e14LOF"
)
n_sample <- length(sample_names)

data_list <- list()
for(i in 1:n_sample){
  data_list[[i]] <- read10X(
    matFile = paste("./data/", sample_names[i], "_matrix.mtx.gz", sep = ""), 
    cellFile = paste("./data/", sample_names[i], "_barcodes.tsv.gz", sep = ""),
    geneFile = "./data/GSE122043_genes.tsv.gz", 
    suffix = sample_names[i], sep = "-"
  )
}
names(data_list) <- sample_names # sparse matrix



# Seurat object for each sample and merge together
GSE122043_mouse_skin <- 
  lapply(data_list, function(x){
  x_S <- CreateSeuratObject(counts =  x, min.features = 1000, names.delim = "-", names.field = 3)
  return(x_S)
}) %>%
  reduce(bind_rows) %>% 
  mutate(sample = orig.ident) %>% 
  mutate(dataset_id = "GSE122043") 
 ## cell type | cell_cluster

GSE122043_mouse_skin %>% saveRDS(GSE122043_mouse_skin,"data/GSE122043_mouse_skin.rds")

data_list[[1]] %>% CreateSeuratObject()
