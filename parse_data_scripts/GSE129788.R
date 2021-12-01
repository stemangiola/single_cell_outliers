### DA-seq on aging mouse brain data
### Original paper: https://www.nature.com/articles/s41593-019-0491-3
### This script reproduces analysis presented in Supp Figure S5,S6
setwd("/stornext/Home/data/allstaff/m/ma.m/single_cell_outliers")
library(Seurat) #V3
library(DAseq)
library(Matrix)
library(reshape2)
library(ggplot2)
library(cowplot)

#source("convenience.R")

## Set Python and GPU
#python2use <- "/data/henry/henry_env/venv/bin/python"
#GPU <- 2

## Set path for FIt-SNE R wrapper
#fitsneR <- "~/git/FIt-SNE/fast_tsne.R"


##=============================================##
## Data prep

## Download data

#if(!dir.exists("./data/")){
#  dir.create("./data/")
#}

# Please go to https://singlecell.broadinstitute.org/single_cell/study/SCP263/aging-mouse-brain#/ to download
# processed data.

# Make sure downloaded data is stored under ./data/ directory.

# Installation (takes time, but only run once):
#install.packages("devtools")
#devtools::install_github("BaderLab/AgingMouseBrainCCInx")

# View the data:
#library(AgingMouseBrainCCInx)
#viewAgingMouseBrainCCInx()

## Load data

# exp mat
data_exp <- read.table(
  "./data/expression_Aging_mouse_brain_portal_data_updated.txt", sep = "\t", header = T, 
  row.names = 1, stringsAsFactors = F
)

data_exp <- Matrix(as.matrix(data_exp), sparse = T)


# meta data
data_meta <- read.table(
  "./data/meta_Aging_mouse_brain_portal_data.txt", sep = "\t", header = T, stringsAsFactors = F
)
data_meta <- data_meta[-1,]
data_meta$cell_type <- gsub("NEUR_mature","mNEUR",data_meta$cell_type)
data_meta$cell_type <- gsub("NEUR_immature","ImmN",data_meta$cell_type)

table(data_meta$cell_type)
table(data_meta$all_cells_by_age)
rownames(data_meta) <- data_meta[,1]

celltype2label <- c(
  "1-OPC","2-OLG","3-OEG","4-NSC","5-ARP","6-ASC","7-NRP","8-ImmN","9-mNEUR","10-NendC",
  "11-EPC","12-HypEPC","13-TNC","14-CPC","15-EC","16-PC","17-VSMC","18-Hb_VC","19-VLMC","20-ABC",
  "21-MG","22-MNC","23-MAC","24-DC","25-NEUT"
)
names(celltype2label) <- sapply(celltype2label, function(x) unlist(strsplit(x,"-"))[2])



## Seurat

# create object
data_S <- CreateSeuratObject(
  counts = data_exp, names.field = 6, names.delim = "_"
)
table(data_S@meta.data$orig.ident)


# add metadata
data_S@meta.data$cell_type <- data_meta[colnames(data_S),"cell_type"]
data_S@meta.data$age <- data_meta[colnames(data_S),"all_cells_by_age"]
data_S@meta.data$condition <- data_S@meta.data$age
data_S@meta.data$condition[data_S@meta.data$age == "2-3mo"] <- "young"
data_S@meta.data$condition[data_S@meta.data$age == "21-22mo"] <- "old"
table(data_S@meta.data$condition)

data_S@meta.data$cell_type_label <- factor(celltype2label[data_S@meta.data$cell_type], levels = celltype2label)
data_S@meta.data$cell_type_num <- as.numeric(data_S@meta.data$cell_type_label)

saveRDS(data_S, "GSE129788_aging_mouse_brain/data/seurat_GSE129788.rds")
# analysis
#data_S <- ScaleData(data_S)
#data_S <- FindVariableFeatures(
#  data_S, selection.method = "mvp", mean.cutoff = c(0.0125,3), dispersion.cutoff = c(0.5,Inf)
#)
#data_S <- RunPCA(
#  data_S, npcs = 20, verbose = F
#)
#data_S <- runFItSNE(
#  data_S, dims.use = 1:20, seed.use = 3, fast.R.path = fitsneR, 
#  ann_not_vptree = FALSE, nthreads = 12
#)

#TSNEPlot(data_S, group.by = "age", pt.size = 0.1)
#TSNEPlot(data_S, group.by = "cell_type_label", label = T, pt.size = 0.1)
