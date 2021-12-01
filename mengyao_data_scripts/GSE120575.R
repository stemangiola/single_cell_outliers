# Defining T cell states associated with response to checkpoint immunotherapy in melanoma
# # https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575 

library(Seurat) #V3
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE120575 
# devtools::install_github("KlugerLab/DAseq")
library(DAseq)
library(Matrix)
library(reshape2)
library(ggplot2)
library(cowplot)
library(tidyseurat)
# install.packages("glue")
# library(glue)
# HCP_scatch = "/stornext/HPCScratch/home/ma.m/single_cell_database/"

# readRDS(glue("{HCP_scatch}/GSE120575_melanoma/GSE120575.R"))

main_S <- readRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/data/seurat_GSE120575.rds")
# source("GSE120575_melanoma_data/convenience.R")

## Set Python and GPU
# python2use <- "/data/henry/henry_env/venv/bin/python"
# GPU <- 3

## Set path for FIt-SNE R wrapper
# fitsneR <- "~/git/FIt-SNE/fast_tsne.R"

##=============================================##
## Data prep

## Load data
# 
# if(!dir.exists("./data/")){
#   dir.create("./data/")
# }

# Expression matrix
# download.file(
#   "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
#   "GSE120575_melanoma_data//data/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
# )

data_exp <- read.table(
  "GSE120575_melanoma/data/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
  sep = "\t", header = F, row.names = 1, stringsAsFactors = F, skip = 2
) 
# header = F:not contains the names of the variable at the first line
# row.names = 1: means the first column gives the row name. 
# character vectors won't be converted to factors.


data_exp <- data_exp[,-16292] # keep all the rows and column without column 16292
data_exp

data_colInfo <- read.table(
  "GSE120575_melanoma/data/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
  sep = "\t", header = F, stringsAsFactors = F, nrows = 2
)

data_colInfo <- data_colInfo[,-1] # keep all the rows and column without column 1


colnames(data_exp) <- data_colInfo[1,] # set the first column into column names
data_exp <- Matrix(as.matrix(data_exp), sparse = T) 


# Patient info
# download.file(
#   "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_patient_ID_single_cells.txt.gz",
#   "GSE120575_melanoma_data/data/GSE120575_patient_ID_single_cells.txt.gz"
# )
patient_info <- read.table(
  "GSE120575_melanoma/data/GSE120575_patient_ID_single_cells.txt.gz", 
  sep = "\t", header = T, stringsAsFactors = F, skip = 19, nrows = 16291
) # nrows: maximum number of rows to read in

patient_info$title
mean(colnames(data_exp) == patient_info$title)
rownames(patient_info) <- patient_info$title # set into row names


# Cluster info
#download.file(
#  "https://raw.githubusercontent.com/KlugerLab/DAseq-paper/master/data/melanoma_cluster_info",
#  "GSE120575_melanoma/data/melanoma_cluster_info"
#)
cluster_info <- read.table(
  "GSE120575_melanoma/data/melanoma_cluster_info", sep = "\t", header = T, stringsAsFactors = F
)
cluster_info$Cell.Name
rownames(cluster_info) <- cluster_info$Cell.Name



## Seurat

data_S <- CreateSeuratObject(
  counts = data_exp, project = "melanoma.immune" # project name
)

colnames(data_S)
# set metadata for each cell

data_S@meta.data$condition <- patient_info[colnames(data_S), "characteristics..response"]
data_S@meta.data$condition # store the response information of the patience info. 
#Ex. Non-responder

data_S@meta.data$lesion <- patient_info[colnames(data_S), 
                                        "characteristics..patinet.ID..Pre.baseline..Post..on.treatment."]
data_S@meta.data$lesion # example output:"Post_P1_2"

data_S@meta.data$cluster <- paste0("G", cluster_info$Cluster.number)
data_S@meta.data$cluster #G1:concatenate vectors

data_S@meta.data$cluster <- factor(data_S@meta.data$cluster, levels = paste0("G", c(1:11)))
# factorize out

cluster2celltype <- c(
  "G1"="G1-B cells", "G2"="G2-Plasma cells", "G3"="G3-Monocytes/Macrophages", "G4"="G4-Dendritic cells",
  "G5"="G5-Lymphocytes", "G6"="G6-Exhausted CD8+ T cells", "G7"="G7-Regulatory T cells", 
  "G8"="G8-Cytotoxicity Lymphocytes", "G9"="G9-Exhausted/HS CD8+ T cells", "G10"="G10-Memory T cells",
  "G11"="G11-Lymphocytes exhausted/cell cycle"
)
cluster2celltype # set the cell type and store

data_S@meta.data$cell_type <- cluster2celltype[as.character(data_S@meta.data$cluster)]
# convert into character data type.
data_S@meta.data$cell_type # Ex.  "G5-Lymphocytes"

data_S@meta.data$cell_type <- factor(data_S@meta.data$cell_type, levels = c(
  cluster2celltype
))

data_S@meta.data$time <- sapply(data_S@meta.data$lesion, FUN = function(x){
  unlist(strsplit(x, split = "_", fixed = T))[1] # split out the characters and only retain the first string character
})
data_S@meta.data$time # Ex.Pre/Post
data_S@meta.data$res_time <- paste(data_S@meta.data$condition, data_S@meta.data$time, sep = "_")
data_S@meta.data$res_time # concatenate together of these two.
# Ex."Non-responder_Pre"
data_S
head(data_S@meta.data)

distinct(data_S,lesion). ## How many samples 

data_S %>% 
  rename(cell_type = cell_type, cell_cluster = cluster,sample = lesion) %>%
  mutate(dataset_id = "GSE120575") %>% 
  saveRDS("/stornext/HPCScratch/home/ma.m/single_cell_database/GSE120575_melanoma/data/GSE120575_melanoma.rds")
