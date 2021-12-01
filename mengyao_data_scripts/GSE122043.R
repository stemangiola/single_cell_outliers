### DA-seq on mouse skin data
### Original paper: https://www.sciencedirect.com/science/article/pii/S1534580718309882
### This script reproduces analysis presented in Figure 3
### https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE122043

# setwd("/stornext/Home/data/allstaff/m/ma.m/single_cell_outliers")
#install.packages("tidyverse")
library(tidyverse)
library(Seurat) #V3
library(DAseq)
library(Matrix)
library(reshape2)
library(ggplot2)
library(cowplot)
library(tidyseurat)

setwd("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE122043_mouse_skin")
# source('convenience.R')

## Set Python and GPU
# python2use <- "/data/henry/henry_env/venv/bin/python"
# GPU <- 2

## Set path for FIt-SNE R wrapper
# fitsneR <- "~/git/FIt-SNE/fast_tsne.R"


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

#if(!dir.exists("./data/")){
#  dir.create("./data/")
#}


# Gene expression files
# download.file(
#  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE122043&format=file", 
#  destfile = "./data/GSE122043_RAW.tar"
#)
system(
  "tar -xvf ./data/GSE122043_RAW.tar -C ./data/"
)

#download.file(
#  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122043/suppl/GSE122043_genes.tsv.gz",
#  "./data/GSE122043_genes.tsv.gz"
#)


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
names(data_list) <- sample_names


# Seurat object for each sample
GSE122043_mouse_skin <- 
  lapply(data_list, function(x){
  x_S <- CreateSeuratObject(counts =  x, min.features = 1000, names.delim = "-", names.field = 3)
  return(x_S)
}) %>%
  reduce(bind_rows) %>% 
  mutate(sample = orig.ident) %>% 
  mutate(dataset_id = "GSE122043") 
 ## cell type | cell_cluster

GSE122043_mouse_skin

saveRDS(GSE122043_mouse_skin,"data/GSE122043_mouse_skin.rds")

'''
  x_S <- subset(x_S, subset = nCount_RNA > 2500 & nCount_RNA < 50000)
  x_S <- NormalizeData(x_S)
  x_S <- ScaleData(x_S)
  x_S <- FindVariableFeatures(x_S, selection.method = "mvp", do.plot = F)
  x_S <- RunPCA(
    x_S, features = rownames(subset(x_S@assays$RNA@meta.features, mvp.dispersion > 0.8)),
    npcs = 10, verbose = F
  )
  x_S <- RunTSNE(x_S, dims = 1:10)
  x_S <- FindNeighbors(x_S, dims = 1:10, verbose = F)
  x_S <- FindClusters(x_S, resolution = 0.1, verbose = F)
  
  return(x_S)
})
# names(data_S_list) <- sample_names
# distinct(data_S_list, orig.ident)

data_S_list %>% 
  counts %>% left_join(annotation_tabe, by="cell") %>% 
  mutate(sample = orig.ident, cell_cluster = seurat_cluster)



lapply(data_S_list, function(x){
  table(x@meta.data$orig.ident)
})
lapply(data_S_list, function(x) DotPlot(x, features = c("Col1a1","Krt14","Krt10"), cols = c("gray","red")))

# select only dermal cells based on Col1a1
dermal_cluster <- lapply(data_S_list, function(x){
  gene.ratio.cluster <- by(x@assays$RNA@data["Col1a1",], INDICES = x@active.ident, 
                           FUN = function(xx) mean(xx > 0))
  gene.ratio.cluster.2 <- by(x@assays$RNA@data["Krt14",] + x@assays$RNA@data["Krt10",], 
                             INDICES = x@active.ident, 
                             FUN = function(xx) mean(xx > 0))
  return(names(gene.ratio.cluster)[gene.ratio.cluster > 0.8 & gene.ratio.cluster.2 < 0.5])
})
data_derm_S_list <- list()
for(i in 1:n_sample){
  data_derm_S_list[[i]] <- subset(
    data_S_list[[i]], cells = which(data_S_list[[i]]@active.ident %in% dermal_cluster[[i]])
  )
}
names(data_derm_S_list) <- sample_names



## Merge data

data_anchors <- FindIntegrationAnchors(object.list = data_derm_S_list)
data_S <- IntegrateData(data_anchors, normalization.method = "LogNormalize")
data_S <- ScaleData(data_S)
data_S <- RunPCA(data_S, verbose = F)
plot(data_S@reductions$pca@stdev)

data_S@meta.data$time <- gsub("Control","",sapply(data_S@meta.data$orig.ident, FUN = function(x){
  unlist(strsplit(x, split = "_", fixed = T))[2]
}))
data_S@meta.data$time <- factor(data_S@meta.data$time, levels = c("e14","e13"))

data_S <- runFItSNE(
  data_S, dims.use = 1:40, seed.use = 3, fast.R.path = fitsneR, 
  ann_not_vptree = FALSE, nthreads = 12
)
TSNEPlot(data_S, group.by = "time")
FeaturePlot(data_S, features = paste("rna",c("Col1a1","Krt14","Krt10"),sep="_"), cols = c("gray","red"))
FeaturePlot(data_S, features = paste("rna",c("Dkk1","Lef1","Ptch1","Bmp4"),sep="_"), cols = c("gray","red"))
''' 
